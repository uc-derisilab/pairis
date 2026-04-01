// Analysis processes

process GROUP_STRUCTURES_BY_BCR {
    label 'grouping'

    input:
    path complexes_dir

    output:
    path "groups.txt"

    script:
    def window_size = params.window_sizes ? params.window_sizes[0] : 15
    def bcr_dir_absolute = file("${workflow.launchDir}/${params.bcr_dir}").toString()
    """
    python3 ${projectDir}/bin/group_structures_by_bcr.py \\
        ${complexes_dir} \\
        ${bcr_dir_absolute} \\
        ${window_size}
    """
}

process RUN_ROSETTA {
    tag "${bcr_name}_${peptide_name}"
    label 'rosetta'

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.af3_output_dir}/results/rosetta/${kmer_size}/${bcr_name}", mode: 'copy'

    input:
    tuple val(kmer_size), val(bcr_name), val(peptide_name), path(complexes_dir)

    output:
    path "${peptide_name}_rosetta.csv"

    script:
    // Pattern must include BCR to filter to only this BCR's structures
    def full_pattern = "${peptide_name}_${bcr_name}"
    """
    python ${projectDir}/bin/rosetta_energy_analysis.py \\
        --input_dir ${complexes_dir} \\
        --pattern "${full_pattern}" \\
        --kmer-size ${kmer_size} \\
        --interface peptide_antibody \\
        --output ${peptide_name}_rosetta.csv
    """
}

process COLLATE_RESULTS {
    publishDir "${params.af3_output_dir}/results/collated", mode: 'copy'

    input:
    val ready

    output:
    path "*.csv"
    path "summary.txt"

    script:
    def af3_input_dir = file("${workflow.launchDir}/${params.outdir}/af3_inputs/complexes").toString()
    """
    python3 ${projectDir}/bin/collate_af3_rosetta_results.py \\
        --af3-input-dir ${af3_input_dir} \\
        --af3-output-dir ${params.af3_output_dir}/af3_outputs/complexes \\
        --rosetta-dir ${params.af3_output_dir}/results/rosetta \\
        --output-dir .

    # Generate summary
    echo "PAIRIS Collated Results Summary" > summary.txt
    echo "===============================" >> summary.txt
    echo "" >> summary.txt
    echo "Output files:" >> summary.txt
    for csv in *.csv; do
        if [ -f "\$csv" ]; then
            rows=\$(tail -n +2 "\$csv" | wc -l)
            echo "  - \$csv: \$rows structures" >> summary.txt
        fi
    done
    """
}
