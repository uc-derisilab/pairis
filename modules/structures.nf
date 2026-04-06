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
    path msa_dir    // staged MSA directory — guarantees files exist before os.path.exists() checks

    output:
    path "*.json"

    script:
    def sliding_args = (window_size && window_size > 0) ? "--sliding -k ${window_size}" : ""
    def msa_arg = msa_index.name != 'NO_FILE' ? "--msa-index ${msa_index}" : ""
    def msa_dir_arg = msa_dir.name != 'NO_MSA_DIR' ? "--msa-dir ${msa_dir}" : ""
    """
    python ${projectDir}/bin/generate_complex_inputs.py \\
        --peptide-fasta ${peptide_fasta} \\
        --antibody-fasta ${bcr} \\
        --antibody-name ${bcr.baseName} \\
        --output-dir . \\
        --num-seeds ${num_seeds} \\
        ${sliding_args} \\
        ${msa_arg} \\
        ${msa_dir_arg}
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
    def is_local = workflow.profile?.contains('local')

    if (is_local) {
        """
        GPU_ID=-1
        for gpu_id in \$(seq 0 \$(( ${params.local_num_gpus} - 1 ))); do
            lockdir="/tmp/pairis_gpu_\${gpu_id}.lock"
            if mkdir "\$lockdir" 2>/dev/null; then
                GPU_ID=\$gpu_id
                trap "rmdir '\$lockdir' 2>/dev/null || true" EXIT
                break
            fi
        done
        if [ \$GPU_ID -eq -1 ]; then
            echo "ERROR: no free GPU slot (all ${params.local_num_gpus} locked)" >&2; exit 1
        fi
        export CUDA_VISIBLE_DEVICES=\$GPU_ID
        mkdir -p ${json.baseName}

        apptainer exec \\
            --env XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \\
            --nv \\
            --bind \${PWD}:/root/input \\
            --bind \${PWD}/${json.baseName}:/root/output \\
            --bind ${params.af3_model_dir}:/root/models \\
            --bind ${params.af3_db_dir}:/root/public_databases \\
            ${params.af3_sif} \\
            python /app/alphafold/run_alphafold.py \\
                --json_path=/root/input/${json} \\
                --model_dir=/root/models \\
                --db_dir=/root/public_databases \\
                --output_dir=/root/output \\
                --flash_attention_implementation=xla \\
                --run_data_pipeline=false

        # Cleanup
        find ${json.baseName} -type d -name "seed*" -exec rm -rf {} + || true
        find ${json.baseName} -type f -name "*_data.json" -exec rm -f {} + || true
        find ${json.baseName} -type f -name "TERMS_OF_USE.md" -exec rm -f {} + || true
        """
    } else {
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
}
