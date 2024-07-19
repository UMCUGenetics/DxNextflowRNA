#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fastq_to_bam (sub)workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow fastq_to_bam {
    take:
        ch_fastq  // channel: [ val(meta), path(bam), path(bai/csi) ]
        ch_star_index  // channel: [ val(meta), path(index)]
        ch_gtf  // channel: [ val(meta), path(gtf)]

    main:
        // Create empty versions channel, and fill with each tools version
        ch_versions = Channel.empty()

        TRIMGALORE(ch_fastq)
        ch_versions.mix(TRIMGALORE.out.versions.first())

        STAR_ALIGN(TRIMGALORE.out.reads, ch_star_index, ch_gtf, false, params.seq_platform, params.seq_center)
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

        SAMTOOLS_MERGE(
            STAR_ALIGN.out.bam_sorted.map {
                meta, bam ->
                    new_id = meta.id.split('_')[0]
                    [ meta + [id: new_id], bam ]
            }.groupTuple(),
            [[ id:'null' ], []],
            [[ id:'null' ], []],
        )
        ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    emit:
        TRIMGALORE.out.reads, emit: trim_reads  // channel: [ val(meta), path(fq.gz) ]
        TRIMGALORE.out.unpaired, emit: trim_unpaired  // channel: [ val(meta), path(fq.gz) ]
        TRIMGALORE.out.html, emit: trim_html  // channel: [ val(meta), path(html) ]
        TRIMGALORE.out.zip, emit: trim_zip  // channel: [ val(meta), path(zip) ]
        TRIMGALORE.out.log, emit: trim_log  // channel: [ val(meta), path(txt) ]

        STAR_ALIGN.out.bam, emit: orig_bam	// channel: [ val(meta), path(bam) ]
        STAR_ALIGN.out.bam_sorted, emit: orig_bam_sorted	// channel: [ val(meta), path(bam) ]
        STAR_ALIGN.out.bam_unsorted, emit: orig_bam_unsorted	// channel: [ val(meta), path(bam) ]
        STAR_ALIGN.out.bam_transcript, emit: orig_bam_transcript	// channel: [ val(meta), path(bam) ]
        STAR_ALIGN.out.sam, emit: sam  // channel: [ val(meta), path(sam) ]
        STAR_ALIGN.out.fastq, emit: fastq	// channel: [ val(meta), path(fastq) ]
        STAR_ALIGN.out.tab, emit: tab	// channel: [ val(meta), path(tab) ]
        STAR_ALIGN.out.spl_junc_tab, emit: spl_junc_tab  // channel: [ val(meta), path(spl_junc_tab) ]
        STAR_ALIGN.out.read_per_gene_tab, emit: read_per_gene_tab  // channel: [ val(meta), path(read_per_gene_tab) ]
        STAR_ALIGN.out.junction, emit: junction  // channel: [ val(meta), path(junction) ]
        STAR_ALIGN.out.log_final, emit: log_final	// channel: [ val(meta), path(log_final) ]
        STAR_ALIGN.out.log_out, emit: log_out	// channel: [ val(meta), path(log_out) ]
        STAR_ALIGN.out.log_progress, emit: log_progress	// channel: [ val(meta), path(log_progress) ]
        STAR_ALIGN.out.wig, emit: wig
        STAR_ALIGN.out.bedgraph, emit: bedgraph

        SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.csi), emit: ch_bam_bai  // channel: [ val(meta), path(bam), path(bai/csi) ]
        SAMTOOLS_MERGE.out.cram.join(SAMTOOLS_MERGE.out.crai), emit: ch_cram_crai  // channel: [ val(meta), path(cram), path(crai) ]

        ch_versions

}