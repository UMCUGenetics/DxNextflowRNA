#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from '../../modules/nf-core/star/align/main'
include { TRIMGALORE } from '../../modules/nf-core/trimgalore/main'

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
        star_ignore_sjdbgtf  // boolean
        seq_platform  // val(seq_center)
        seq_center  // val(seq_center)

    main:
        // Create empty versions channel, and fill with each tools version
        versions = Channel.empty()

        TRIMGALORE(ch_fastq)
        versions.mix(TRIMGALORE.out.versions.first())

        STAR_ALIGN(TRIMGALORE.out.reads, ch_star_index, ch_gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
        versions = versions.mix(STAR_ALIGN.out.versions.first())

        SAMTOOLS_MERGE(
            STAR_ALIGN.out.bam_sorted.map {
                meta, bam ->
                    new_id = meta.id.split('_')[0]
                    [ meta + [id: new_id], bam ]
            }.groupTuple(),
            [[ id:'null' ], []],
            [[ id:'null' ], []],
        )
        versions.mix(SAMTOOLS_MERGE.out.versions.first())

    emit:
        trim_reads= TRIMGALORE.out.reads  // channel: [ val(meta), path(fq.gz) ]
        trim_unpaired= TRIMGALORE.out.unpaired  // channel: [ val(meta), path(fq.gz) ]
        trim_html= TRIMGALORE.out.html  // channel: [ val(meta), path(html) ]
        trim_zip= TRIMGALORE.out.zip  // channel: [ val(meta), path(zip) ]
        trim_log= TRIMGALORE.out.log  // channel: [ val(meta), path(txt) ]

        star_align_bam= STAR_ALIGN.out.bam  // channel: [ val(meta), path(bam) ]
        star_align_bam_sorted= STAR_ALIGN.out.bam_sorted  // channel: [ val(meta), path(bam) ]
        star_align_bam_unsorted= STAR_ALIGN.out.bam_unsorted  // channel: [ val(meta), path(bam) ]
        star_align_bam_transcript= STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), path(bam) ]
        star_align_sam = STAR_ALIGN.out.sam  // channel: [ val(meta), path(sam) ]
        star_align_fastq = STAR_ALIGN.out.fastq  // channel: [ val(meta), path(fastq) ]
        star_align_tab = STAR_ALIGN.out.tab  // channel: [ val(meta), path(tab) ]
        star_align_spl_junc_tab = STAR_ALIGN.out.spl_junc_tab  // channel: [ val(meta), path(spl_junc_tab) ]
        star_align_read_per_gene_tab = STAR_ALIGN.out.read_per_gene_tab  // channel: [ val(meta), path(read_per_gene_tab) ]
        star_align_junction = STAR_ALIGN.out.junction  // channel: [ val(meta), path(junction) ]
        star_align_log_final = STAR_ALIGN.out.log_final  // channel: [ val(meta), path(log_final) ]
        star_align_log_out = STAR_ALIGN.out.log_out  // channel: [ val(meta), path(log_out) ]
        star_align_log_progress = STAR_ALIGN.out.log_progress  // channel: [ val(meta), path(log_progress) ]
        star_align_wig = STAR_ALIGN.out.wig
        star_align_bedgraph = STAR_ALIGN.out.bedgraph

        ch_bam_bai = SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_MERGE.out.csi)  // channel: [ val(meta), path(bam), path(bai/csi) ]
        ch_cram_crai = SAMTOOLS_MERGE.out.cram.join(SAMTOOLS_MERGE.out.crai)  // channel: [ val(meta), path(cram), path(crai) ]

        versions

}
