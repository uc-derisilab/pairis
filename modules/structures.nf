// Structure prediction processes

process GENERATE_COMPLEX_INPUTS {
    tag "${bcr.baseName}"
    label 'complex_input'

    publishDir "${params.outdir}/af3_inputs/complexes", mode: 'copy'

    input:
    path peptide_fasta
    path bcr
    val window_size  // 0 for full sequence, or positive integer for sliding windows
    val num_seeds
    path msa_index  // MSA index file (optional)

    output:
    path "*.json"

    script:
    def sliding_args = (window_size && window_size > 0) ? "--sliding -k ${window_size}" : ""
    def msa_arg = msa_index.name != 'NO_FILE' ? "--msa-index ${msa_index}" : ""
    """
    python ${projectDir}/bin/generate_complex_inputs.py \\
        --peptide-fasta ${peptide_fasta} \\
        --antibody-fasta ${bcr} \\
        --antibody-name ${bcr.baseName} \\
        --output-dir . \\
        --num-seeds ${num_seeds} \\
        ${sliding_args} \\
        ${msa_arg}
    """
}

process RUN_AF3_FOLDING {
    tag "${json.baseName}"
    label 'af3_folding'

    errorStrategy 'retry'
    maxRetries 10

    publishDir "${params.af3_output_dir}/af3_outputs/complexes", mode: 'copy'

    input:
    path json

    output:
    path "${json.baseName}/*/*"

    script:
    """
    module load ${params.af3_module}

    mkdir -p ${json.baseName}

    alphafold \\
        --json_path=${json} \\
        --output_dir=${json.baseName} \\
        --run_data_pipeline=false

    # Remove seed folders from AF3 run to save space
    find ${json.baseName} -type d -name "seed*" -exec rm -rf {} +

    # Remove _data.json files to save space (already saved from precomputed MSAs)
    find ${json.baseName} -type f -name "*_data.json" -exec rm -f {} +

    # Remove TERMS_OF_USE.md files
    find ${json.baseName} -type f -name "TERMS_OF_USE.md" -exec rm -f {} +
    """
}
