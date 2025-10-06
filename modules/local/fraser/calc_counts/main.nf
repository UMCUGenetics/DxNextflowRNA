process FRASER_CALC_COUNTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.ag.eg.db_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:3d3e34783fd3bcc9' :
        'community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:6ff8a97056de00fa'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_junction_counts.tsv.gz"), emit: junctionCounts

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ""
    """

    get_fraser_counts.R \\
        --refset ${bam} \\
        --threads ${task.cpus} \\
        --outdir ./
    """
}
