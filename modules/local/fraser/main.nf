
process FRASER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "oras://community.wave.seqera.io/library/bioconductor-fraser_r-argparse:f6da91bb26f0dbed"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(ref_bam), path(ref_bai)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*_heatmap.pdf"), emit: heatmap
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fraser.R \\
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
}
