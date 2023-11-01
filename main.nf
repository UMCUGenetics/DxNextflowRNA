#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowRNA
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from './modules/nf-core/fastqc/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'

include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    ch_fastq = Channel.fromFilePairs("$params.input/*_R{1,2}_001.fastq.gz")
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    //TRIMGALORE

    FASTQC(ch_fastq)

    SAMTOOLS_INDEX ( STAR.out.bam_sorted )
    STAR.out.bam_sorted
        .join(SAMTOOLS_INDEX.out.bai)
        .set { ch_bam_bai }

    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS( ch_bam_bai )

    SUBREAD_FEATURECOUNTS(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam.map{ meta, bam -> [ meta, bam, params.genome ]})

}
