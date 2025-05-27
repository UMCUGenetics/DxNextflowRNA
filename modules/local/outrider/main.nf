process OUTRIDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/bioconductor-outrider_r-argparse_r-tidyverse:bf4b6d86b39e3c3e"


    input:
    tuple path(counts), path(refset)
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.outrider_result.tsv"), emit: tsv
    tuple val(meta), path("counts_heatplots.pdf"), emit: psv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def nthreads = task.cpus ?: 1
    def count_files = counts.join(" ")
    """
    outrider.R \\
        ${args} \\
        --ref ${refset} \\
        --pref ${prefix} \\
        --gtf ${gtf} \\
        --threads ${nthreads} \\
        ${count_files} \\
    """
}
