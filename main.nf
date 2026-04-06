#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { GENERATE_MSA_INPUTS as GEN_MSA_PEPTIDES } from './modules/msa'
include { GENERATE_MSA_INPUTS as GEN_MSA_BCRS } from './modules/msa'
include { RUN_AF3_MSA } from './modules/msa'
include { RUN_ALPHAFAST_MSA } from './modules/msa'
include { EXTRACT_MSAS } from './modules/msa'
include { GENERATE_COMPLEX_INPUTS } from './modules/structures'
include { RUN_AF3_FOLDING } from './modules/structures'
include { GROUP_STRUCTURES_BY_BCR } from './modules/analysis'
include { RUN_ROSETTA } from './modules/analysis'
include { COLLATE_RESULTS } from './modules/analysis'

workflow {
    // Validate local profile parameters
    if (workflow.profile?.contains('local')) {
        if (!params.af3_sif) {
            error "ERROR: -profile local requires --af3_sif (path to alphafold.sif)"
        }
        if (!params.af3_model_dir) {
            error "ERROR: -profile local requires --af3_model_dir (path to AF3 model params)"
        }
        if (!params.af3_db_dir) {
            error "ERROR: -profile local requires --af3_db_dir (path to AF3 databases)"
        }
        if (!params.local_venv) {
            log.warn "WARNING: --local_venv not set. Python scripts may fail if python is not on PATH."
        }
    }

    // Validate AlphaFast params (applies to all profiles)
    if (params.msa_engine == 'alphafast') {
        if (!params.alphafast_sif) {
            error "ERROR: msa_engine='alphafast' requires --alphafast_sif (path to alphafast.sif)"
        }
        if (!params.alphafast_mmseqs_db_dir) {
            error "ERROR: msa_engine='alphafast' requires --alphafast_mmseqs_db_dir (path to MMseqs2 padded databases)"
        }
    }

    // Load input files
    peptide_fasta = file(params.peptide_fasta)
    bcr_files = Channel.fromPath("${params.bcr_dir}/*.fasta")

    // Initialize structures channel
    structures = Channel.empty()

    // ===== Phase 1: MSA Generation (optional) =====
    if (params.run_msa_generation) {
        // Generate MSA inputs for peptides with different window sizes
        if (params.window_sizes) {
            peptide_msa_jsons = GEN_MSA_PEPTIDES(
                Channel.value(peptide_fasta),
                Channel.from(params.window_sizes)
            )
        } else {
            // No sliding window, use full sequences
            peptide_msa_jsons = GEN_MSA_PEPTIDES(peptide_fasta, 0)
        }

        // Generate MSA inputs for BCRs (no sliding windows, use 0 to indicate full sequence)
        // Use chain_dir (individual chain FASTAs) if available, otherwise fall back to bcr_dir
        def msa_input_dir = params.chain_dir ?: params.bcr_dir
        chain_files = Channel.fromPath("${msa_input_dir}/*.fasta")
        bcr_msa_jsons = GEN_MSA_BCRS(chain_files, 0)

        // Combine all MSA inputs
        all_msa_jsons = peptide_msa_jsons
            .flatten()
            .mix(bcr_msa_jsons.flatten())

        // Run MSA generation (engine selected by params.msa_engine)
        if (params.msa_engine == 'alphafast') {
            msa_outputs = RUN_ALPHAFAST_MSA(all_msa_jsons)
        } else {
            msa_outputs = RUN_AF3_MSA(all_msa_jsons)
        }

        // Extract MSAs and build index
        extraction_results = EXTRACT_MSAS(msa_outputs.collect())
        msa_index_file = extraction_results.index
        msa_dir = extraction_results.msa_dir
    } else {
        // Load existing MSA index if available
        def msa_index_path = "${params.af3_output_dir}/results/msas/msa_index.json"
        def index_file = file(msa_index_path).exists() ? file(msa_index_path) : file('NO_FILE')
        msa_index_file = Channel.value(index_file)
        def msa_dir_path = file("${params.af3_output_dir}/results/msas/msas")
        msa_dir = Channel.value(msa_dir_path.exists() ? msa_dir_path : file('NO_MSA_DIR'))
    }

    // ===== Phase 2: Structure Prediction =====
    if (params.run_structure_prediction) {
        // Generate complex JSONs for each BCR with peptides
        if (params.window_sizes) {
            inputs = bcr_files
                .combine(Channel.from(params.window_sizes))
                .combine(msa_index_file)
                .combine(msa_dir)
                .multiMap { bcr, ws, msa_idx, msa_d ->
                    peptide: peptide_fasta
                    bcr: bcr
                    window: ws
                    seeds: params.num_seeds
                    index: msa_idx
                    msa_directory: msa_d
                }

            complex_jsons = GENERATE_COMPLEX_INPUTS(
                inputs.peptide,
                inputs.bcr,
                inputs.window,
                inputs.seeds,
                inputs.index,
                inputs.msa_directory
            )
        } else {
            // No sliding window, use full sequences
            complex_jsons = GENERATE_COMPLEX_INPUTS(
                peptide_fasta,
                bcr_files,
                0,
                params.num_seeds,
                msa_index_file,
                msa_dir
            )
        }

        // Run AF3 for structure prediction
        structures = RUN_AF3_FOLDING(complex_jsons.flatten())
    } else {
        // Load existing structures from previous run
        structures = Channel.fromPath("${params.af3_output_dir}/af3_outputs/complexes/*/*/*_model*.cif")
    }

    // ===== Phase 3: Analysis =====
    if (params.run_rosetta_analysis) {
        // Wait for all structures to complete, then pass directory for grouping
        complexes_dir = structures.collect()
            .map { file("${params.af3_output_dir}/af3_outputs/complexes") }

        // Group structures by (kmer_size, BCR, peptide)
        groups_file = GROUP_STRUCTURES_BY_BCR(complexes_dir)

        // Parse groups file and create channel (kmer_size, BCR, peptide)
        groups_info = groups_file
            .splitText()
            .map { line ->
                def parts = line.trim().split(',')
                tuple(parts[0], parts[1], parts[2])
            }

        // Combine with complexes_dir to create full tuple
        bcr_groups = groups_info.combine(complexes_dir)

        // Process each BCR-peptide group
        rosetta_results = RUN_ROSETTA(bcr_groups)

        // Collate AF3 iPTM scores and Rosetta binding energies
        collated_results = COLLATE_RESULTS(rosetta_results.collect().map { true })
    }
}

workflow.onComplete {
    println ""
    println "=========================================="
    println "PAIRIS Pipeline Complete!"
    println "=========================================="
    println "Results directory: ${params.outdir}"
    if (params.run_msa_generation) {
        println "AF3 MSA outputs: ${params.af3_output_dir}/af3_outputs/msa_only"
    }
    if (params.run_structure_prediction) {
        println "AF3 structures: ${params.af3_output_dir}/af3_outputs/complexes"
    }
    println ""

    // Generate timing summary from trace.txt.
    // Sleep briefly so the async trace writer has time to flush the last records.
    try {
        Thread.sleep(2000)
        def traceFile   = "${params.outdir}/reports/trace.txt"
        def summaryFile = "${params.outdir}/reports/timing_summary.txt"
        def cmd = [
            "python3",
            "${projectDir}/bin/generate_timing_summary.py",
            "--trace", traceFile,
            "--output", summaryFile,
            "--pipeline-duration", workflow.duration.toString()
        ]
        def proc = cmd.execute()
        proc.waitFor()
        def out = proc.text.trim()
        if (out) {
            println out
            println ""
        }
        println "Timing summary written to: ${summaryFile}"
        println ""
    } catch (Exception e) {
        log.warn "Could not generate timing summary: ${e.message}"
    }
}
