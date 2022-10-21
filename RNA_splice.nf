#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq modules
//include extractFastqPairFromDir from './modules/fastq.nf'
//def fastq_files = extractFastqPairFromDir(params.fastq_path)

//QC module
//include FastQC from './modules/FastQC.nf' params(optional:'')

//AS calling_rmats
include RMATS from './modules/rmats.nf'

//ASGL ## ME for variantcalling
include ASGAL from './modules/ASGAL.nf' params(optional: "--multi")


transcriptome_ch = Channel.fromPath(params.transcripts)
genome_ch = Channel.fromPath(params.genome)
gtf_ch = Channel.fromPath(params.gtf)
STAR_genome_index_ch = Channel.fromPath(params.STAR_genome_index)
//vcf_file_ch=Channel.fromPath(params.vcf_infile)
//bam_infile_ch=Channel.fromPath(params.bam_infile)

workflow {
    // FASTQC
    //fastqc_ch=FastQC(fastq_files)

    //RMATS
    rmats_AS_calling_ch=RMATS(params.s1, gtf_ch, STAR_genome_index_ch)

    //Index BAM file needed here
 
    //ASGAL
    //ASGAL_ch=ASGAL(genome_ch, transcriptome_ch, gtf_ch, fastq_files)
}
