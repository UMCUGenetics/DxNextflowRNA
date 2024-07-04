#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fastq_to_bam (sub)workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow fastq_to_bam {
    take: ch_fastq

    main:
        // Reference file channels
        ch_star_index = Channel.fromPath(params.star_index).map {star_index -> [star_index.getSimpleName(), star_index] }.first()
        ch_gtf = Channel.fromPath(params.gtf).map { gtf -> [gtf.getSimpleName(), gtf] }.first()

        TRIMGALORE( ch_fastq )

        STAR_ALIGN(
            TRIMGALORE.out.reads,
            ch_star_index,
            ch_gtf,
            false,
	        params.seq_platform,
            params.seq_center
        )

        SAMTOOLS_MERGE(
            STAR_ALIGN.out.bam_sorted.map {
                meta, bam ->
                    new_id = meta.id.split('_')[0]
                    [ meta + [id: new_id], bam ]
            }.groupTuple(),
            [ [ id:'null' ], []],
            [ [ id:'null' ], []],
        )

        SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

        SAMTOOLS_MERGE.out.bam
            .join(SAMTOOLS_INDEX.out.bai)
            .set { ch_bam_bai }

    emit:
	ch_bam_bai
}