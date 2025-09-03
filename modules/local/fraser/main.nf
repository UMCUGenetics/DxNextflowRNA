
process FRASER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.ag.eg.db_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:3d3e34783fd3bcc9' :
        'community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:6ff8a97056de00fa'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(ref_bam), path(ref_bai)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """
    fraser.R \\
        ${args} \\
        --input ${bam} \\
        --refset ${ref_bam} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --paired TRUE

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        fraser: \$(fraser.R --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        fraser: \$(fraser.R --version)
    END_VERSIONS
    """
}
