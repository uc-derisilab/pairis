// MSA generation processes

process GENERATE_MSA_INPUTS {
    tag "${fasta.baseName}"
    label 'msa_input'

    publishDir "${params.outdir}/af3_inputs/msa", mode: 'copy'

    input:
    path fasta
    val window_size  // 0 for full sequence, or positive integer for sliding windows

    output:
    path "*.json"

    script:
    def sliding_args = (window_size && window_size > 0) ? "--sliding -k ${window_size}" : ""
    """
    python ${projectDir}/bin/generate_msa_af3_input.py \\
        --fasta ${fasta} \\
        --output-dir . \\
        ${sliding_args}
    """
}

process RUN_AF3_MSA {
    tag "${json.baseName}"
    label 'af3_msa'

    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.af3_output_dir}/af3_outputs/msa_only", mode: 'copy'

    input:
    path json

    output:
    path "${json.baseName}/*/*_data.json"

    script:
    """
    module load ${params.af3_module}

    mkdir -p ${json.baseName}

    alphafold \\
        --json_path=${json} \\
        --output_dir=${json.baseName} \\
        --run_inference=false
    """
}

process EXTRACT_MSAS {
    label 'extract'

    publishDir "${params.af3_output_dir}/results/msas", mode: 'copy'

    input:
    path data_jsons

    output:
    path "msas/*/*.a3m", emit: msa_files
    path "msas/*/metadata.json", emit: metadata_files
    path "msas/*/templates/*.cif", optional: true, emit: template_files
    path "msa_index.json", emit: index

    script:
    """
    # Extract MSAs to organized directory structure
    python ${projectDir}/bin/extract_msas_from_precomputed.py \\
        --input-dir . \\
        --output-dir msas

    # Build index mapping sequence IDs to MSA paths (using published location)
    python ${projectDir}/bin/build_msa_index.py \\
        --msa-dir msas \\
        --published-base ${params.af3_output_dir}/results/msas/msas \\
        --output msa_index.json
    """
}
