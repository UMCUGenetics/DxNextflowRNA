#!/usr/bin/env nextflow
nextflow.preview.dsl=2


params.reads = '/hpc/diaggen/users/Behzad/bam/PMABM000_fastq/*_r{1,2}.fastq'
params.transcripts = '/hpc/diaggen/users/Behzad/bam/STAR/hg38_genome/gencode.v39.transcripts.fa'


process INDEX {

    input:
    path fasta

    output:
    path 'index'

    """
    salmon index \\
        -t "${fasta}" \\
        -i index \\
    """
}


transcriptome_ch = Channel.fromPath(params.transcripts)


process QUANT {

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path (pair_id)

    """
    salmon quant \\
        --index $index\\
        --libType=U  \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o $pair_id
    """
}


read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

workflow {

    index_ch=INDEX(transcriptome_ch)
    quant_ch=QUANT(index_ch,read_pairs_ch)

}
