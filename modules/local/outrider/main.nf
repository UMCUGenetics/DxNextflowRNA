process OUTRIDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-outrider_r-argparse_r-tidyverse:bf4b6d86b39e3c3e' :
        'community.wave.seqera.io/library/bioconductor-outrider_r-argparse_r-tidyverse:7658fd7a410ca8f5' }"



    input:
    tuple val(meta), path(query), path(refset)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*_outrider_result_full.tsv")      , emit: tsv_full
    tuple val(meta), path("*_outrider_result_signif.tsv")    , emit: tsv_signif
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def nthreads = task.cpus ?: 1
    """
    outrider.R \\
        ${args} \\
        --ref ${refset} \\
        --gtf ${gtf} \\
        --threads ${nthreads} \\
        --queries ${query} \\
        --prefix ${prefix}

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        outrider: \$(outrider.R --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_outrider_result_full.tsv
    touch ${prefix}_outrider_result_signif.tsv

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        outrider: \$(outrider.R --version)
    END_VERSIONS
    """
}
