process FRASER_BUILDCOUNTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.ag.eg.db_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:3d3e34783fd3bcc9' :
        'community.wave.seqera.io/library/bioconductor-fraser_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-argparse:6ff8a97056de00fa'}"

    input:
    tuple val(meta), path(refset_junctions)
    tuple val(meta2), path(sample_junctions)

    output:
    tuple val(meta), path("*_countTable.tsv"),           emit: countTable
    tuple val(meta), path("*_spliceSite_counts.tsv.gz"), emit: spliceCounts
    tuple val(meta), path("*_junction_counts.tsv.gz"),   emit: junctionCounts


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """

    fraser_buildcounts.R \\
        ${args} \\
        --prefix ${prefix} \\
        --refset_junctions ${refset_junctions} \\
        --sample_junctions ${sample_junctions}

    """

}
