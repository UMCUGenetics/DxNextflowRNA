process OUTRIDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/bioconductor-outrider_r-argparse_r-tidyverse:bf4b6d86b39e3c3e"


    input:
    tuple val(meta), path(query), path(refset)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*_outrider_result_full.tsv")      , emit: tsv_full
    tuple val(meta), path("*_outrider_result_signif.tsv")    , emit: tsv_signif
    tuple val(meta), path("*_outrider_result_signif_mqc.tsv"), emit: multiqc_tsv
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


    # Add metadata to the _mqc.tsv so multiqc knows how to format them in the report
    echo "# parent_id: OUTRIDER - ${task.ext.mode}" >> ${prefix}_outrider_result_signif_mqc.tsv
    echo "# parent_name: OUTRIDER - ${task.ext.mode}" >> ${prefix}_outrider_result_signif_mqc.tsv
    echo "# section_name: ${prefix}" >> ${prefix}_outrider_result_signif_mqc.tsv
    cat ${prefix}_outrider_result_signif.tsv >> ${prefix}_outrider_result_signif_mqc.tsv

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        outrider: \$(outrider.R --version)
    END_VERSIONS
    """
}
