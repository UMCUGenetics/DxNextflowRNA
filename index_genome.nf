#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Index genome for STAR
include GenomeGenerate from './modules/GenomeGenerate.nf'

//transcriptome_ch = Channel.fromPath(params.transcripts)
genome_ch = Channel.fromPath(params.genome)
gtf_ch = Channel.fromPath(params.gtf)


workflow {

    //Index genome for STAR    
    GenomeGenerate(genome_ch, gtf_ch)

}
