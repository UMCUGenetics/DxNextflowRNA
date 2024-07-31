#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE } from '../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_EXON } from '../modules/nf-core/subread/featurecounts/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    featurecounts subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow featurecounts {
    take:
        ch_bam_bai    
    main:
        SUBREAD_FEATURECOUNTS_GENE(
            ch_bam_bai.map { meta, bam, bai -> [meta, bam, params.gtf] }
        )

        SUBREAD_FEATURECOUNTS_EXON(
            ch_bam_bai.map { meta, bam, bai -> [meta, bam, params.gtf] }
        )

    emit:
	    SUBREAD_FEATURECOUNTS_GENE.out.counts
	    SUBREAD_FEATURECOUNTS_EXON.out.counts
}


workflow featurecounts_entry {
    ch_bam_bai = Channel.fromFilePairs("$params.input/*.{bam,bam.bai}").map { meta, bambai -> [[id:bambai[0].getSimpleName()], bambai[0], bambai[1]] }
    featurecounts(ch_bam_bai)
}