#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    featurecounts (sub)workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow featurecounts {
    take:
        ch_bam_bai    
    main:
        ch_bam_gene = ch_bam_bai.map { meta, bam, bai -> [meta, bam, params.gtf, "gene"] }
        ch_bam_exon = ch_bam_bai.map { meta, bam, bai -> [meta, bam, params.gtf, "exon"] }

        SUBREAD_FEATURECOUNTS(
            ch_bam_gene.concat(ch_bam_exon)
        )

    emit:
	SUBREAD_FEATURECOUNTS.out.counts
}

workflow featurecounts_entry {
    ch_bam_bai = Channel.fromFilePairs("$params.outdir/bam_files/*.{bam,bam.bai}").map { meta, bambai -> [[id:bambai[0].getSimpleName()], bambai[0], bambai[1]] }
    featurecounts(ch_bam_bai)
}