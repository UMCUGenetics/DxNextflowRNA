#!/usr/bin/env nextflow
nextflow.preview.dsl=2


// Indexing modules
include { INDEX } from './modules/salmon_index.nf'

// Quantification modules
include { QUANT } from './modules/salmon_quant.nf'

// fastq modules
include extractFastqPairFromDir from './modules/fastq.nf'

def fastq_files = extractFastqPairFromDir(params.fastq_path)


transcriptome_ch = Channel.fromPath(params.transcripts)


//read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

workflow {

    index_ch=INDEX(transcriptome_ch)
    quant_ch=QUANT(index_ch,fastq_files)

}
