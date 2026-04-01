#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { GENERATE_MSA_INPUTS as GEN_MSA_PEPTIDES } from './modules/msa'
include { GENERATE_MSA_INPUTS as GEN_MSA_BCRS } from './modules/msa'
include { RUN_AF3_MSA } from './modules/msa'
include { EXTRACT_MSAS } from './modules/msa'
include { GENERATE_COMPLEX_INPUTS } from './modules/structures'
include { RUN_AF3_FOLDING } from './modules/structures'
include { GROUP_STRUCTURES_BY_BCR } from './modules/analysis'
include { RUN_ROSETTA } from './modules/analysis'
include { COLLATE_RESULTS } from './modules/analysis'

workflow {
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

        // Run AF3 for MSA generation
        msa_outputs = RUN_AF3_MSA(all_msa_jsons)

        // Extract MSAs and build index
        extraction_results = EXTRACT_MSAS(msa_outputs.collect())
        msa_index_file = extraction_results.index
    } else {
        // Load existing MSA index if available
        def msa_index_path = "${params.af3_output_dir}/results/msas/msa_index.json"
        def index_file = file(msa_index_path).exists() ? file(msa_index_path) : file('NO_FILE')
        msa_index_file = Channel.value(index_file)
    }

    // ===== Phase 2: Structure Prediction =====
    if (params.run_structure_prediction) {
        // Generate complex JSONs for each BCR with peptides
        if (params.window_sizes) {
            inputs = bcr_files
                .combine(Channel.from(params.window_sizes))
                .combine(msa_index_file)
                .multiMap { bcr, ws, msa_idx ->
                    peptide: peptide_fasta
                    bcr: bcr
                    window: ws
                    seeds: params.num_seeds
                    index: msa_idx
                }

            complex_jsons = GENERATE_COMPLEX_INPUTS(
                inputs.peptide,
                inputs.bcr,
                inputs.window,
                inputs.seeds,
                inputs.index
            )
        } else {
            // No sliding window, use full sequences
            complex_jsons = GENERATE_COMPLEX_INPUTS(
                peptide_fasta,
                bcr_files,
                0,
                params.num_seeds,
                msa_index_file
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
}
