process OUTRIDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/enfironment.yml"
    container

    input:
    tuple val(meta), path(counts)
    tuple val(meta2), path(refset)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.pdf"), emit: psv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def nthreads = task.cpus ?: 1
    """
    outrider.R \\
        ${args} \\
        ${counts} \\
        --ref ${refset} \\
        --pref ${prefix} \\
        --gtf ${gtf} \\
        --threads %{nthreads}
    """
}
