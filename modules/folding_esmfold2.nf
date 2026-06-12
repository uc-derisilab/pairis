// ESMFold2 structure prediction processes
//
// Mirrors modules/structures.nf (the AlphaFold3 backend) so that ESMFold2 is a
// drop-in folding backend: inputs publish to the SAME af3_inputs/complexes dir
// and folded outputs land in the SAME af3_outputs/complexes layout that the
// downstream Rosetta + collation stages consume unchanged.

process GENERATE_ESMFOLD2_INPUTS {
    tag "${bcr.baseName}"
    label 'complex_input'

    publishDir "${params.outdir}/af3_inputs/complexes", mode: 'copy'

    input:
    path peptide_fasta
    path bcr
    val window_size  // 0 for full sequence, or positive integer for sliding windows
    path msa_index   // MSA index file (optional)
    path msa_dir     // staged MSA directory — guarantees files exist before os.path.exists() checks

    output:
    path "*.json"

    script:
    def sliding_args = (window_size && window_size > 0) ? "--sliding -k ${window_size}" : ""
    def msa_arg = msa_index.name != 'NO_FILE' ? "--msa-index ${msa_index}" : ""
    def msa_dir_arg = msa_dir.name != 'NO_MSA_DIR' ? "--msa-dir ${msa_dir}" : ""
    """
    python ${projectDir}/bin/generate_esmfold2_inputs.py \\
        --peptide-fasta ${peptide_fasta} \\
        --antibody-fasta ${bcr} \\
        --antibody-name ${bcr.baseName} \\
        --output-dir . \\
        --model ${params.esmfold2_model} \\
        --num-loops ${params.esmfold2_num_loops} \\
        --num-sampling-steps ${params.esmfold2_num_sampling_steps} \\
        ${sliding_args} \\
        ${msa_arg} \\
        ${msa_dir_arg}
    """
}

process RUN_ESMFOLD2_FOLDING {
    tag "${json.baseName}"
    label 'esmfold2_folding'

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

        python ${projectDir}/bin/run_esmfold2.py \\
            --input ${json} \\
            --output-dir ${json.baseName} \\
            --model ${params.esmfold2_model} \\
            --num-seeds ${params.num_seeds}
        """
    } else {
        """
        mkdir -p ${json.baseName}

        python ${projectDir}/bin/run_esmfold2.py \\
            --input ${json} \\
            --output-dir ${json.baseName} \\
            --model ${params.esmfold2_model} \\
            --num-seeds ${params.num_seeds}
        """
    }
}
