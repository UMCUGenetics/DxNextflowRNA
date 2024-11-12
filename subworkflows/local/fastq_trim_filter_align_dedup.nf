#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_CONVERT                  } from '../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX                    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                    } from '../../modules/nf-core/samtools/merge/main'
include { SORTMERNA as SORTMERNA_READS      } from '../../modules/local/sortmerna/main'

include { STAR_ALIGN                        } from '../../modules/nf-core/star/align/main'
include { TRIMGALORE                        } from '../../modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import subworkflows, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from '../../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ_TRIM_FILTER_ALIGN_DEDUP (sub)workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FASTQ_TRIM_FILTER_ALIGN_DEDUP {
    take:
    ch_fasta_fai        // channel: [ val(meta), path(fa), path(fai) ]
    ch_fastq            // channel: [ val(meta), [ path(fastq1), path(fastq2) ] ]
    ch_gtf              // channel: [ val(meta), path(gtf) ]
    ch_sortmerna_fastas // channel: [ val(meta), file(fastas)]
    ch_sortmerna_index  // channel: [ val(meta), path(/path/to/sortmerna/index/)
    ch_star_index       // channel: [ val(meta), path(star_index) ]
    seq_platform        // val(seq_platform)
    seq_center          // val(seq_center)
    star_ignore_sjdbgtf // boolean

    main:
    // Create empty versions channel, and fill with each tools version
    ch_versions = Channel.empty()

    TRIMGALORE(ch_fastq)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    SORTMERNA_READS(TRIMGALORE.out.reads, ch_sortmerna_fastas, ch_sortmerna_index)
    ch_versions = ch_versions.mix(SORTMERNA_READS.out.versions.first())

    STAR_ALIGN(SORTMERNA_READS.out.reads, ch_star_index, ch_gtf, star_ignore_sjdbgtf, seq_platform, seq_center)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SAMTOOLS_MERGE(
        STAR_ALIGN.out.bam_sorted.map { meta, bam ->
            new_id = meta.id.split('_')[0]
            [meta + [id: new_id], bam]
        }.groupTuple(),
        [[id: 'null'], []],
        [[id: 'null'], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS(SAMTOOLS_MERGE.out.bam.join(SAMTOOLS_INDEX.out.bai), true)
    ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.versions)

    SAMTOOLS_CONVERT(
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam.join(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bai),
        ch_fasta_fai.map { meta, fasta, fai -> [meta, fasta] },
        ch_fasta_fai.map { meta, fasta, fai -> [meta, fai] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    emit:
    trim_reads                   = TRIMGALORE.out.reads // channel: [ val(meta), path(fq.gz) ]
    trim_unpaired                = TRIMGALORE.out.unpaired // channel: [ val(meta), path(fq.gz) ]
    trim_html                    = TRIMGALORE.out.html // channel: [ val(meta), path(html) ]
    trim_zip                     = TRIMGALORE.out.zip // channel: [ val(meta), path(zip) ]
    trim_log                     = TRIMGALORE.out.log // channel: [ val(meta), path(txt) ]

    sortmerna_log                = SORTMERNA_READS.out.log // channel: [ val(meta), path(log) ]

    star_align_bam               = STAR_ALIGN.out.bam // channel: [ val(meta), path(bam) ]
    star_align_bam_sorted        = STAR_ALIGN.out.bam_sorted // channel: [ val(meta), path(bam) ]
    star_align_bam_unsorted      = STAR_ALIGN.out.bam_unsorted // channel: [ val(meta), path(bam) ]
    star_align_bam_transcript    = STAR_ALIGN.out.bam_transcript // channel: [ val(meta), path(bam) ]
    star_align_sam               = STAR_ALIGN.out.sam // channel: [ val(meta), path(sam) ]
    star_align_fastq             = STAR_ALIGN.out.fastq // channel: [ val(meta), path(fastq) ]
    star_align_tab               = STAR_ALIGN.out.tab // channel: [ val(meta), path(tab) ]
    star_align_spl_junc_tab      = STAR_ALIGN.out.spl_junc_tab // channel: [ val(meta), path(spl_junc_tab) ]
    star_align_read_per_gene_tab = STAR_ALIGN.out.read_per_gene_tab // channel: [ val(meta), path(read_per_gene_tab) ]
    star_align_junction          = STAR_ALIGN.out.junction // channel: [ val(meta), path(junction) ]
    star_align_log_final         = STAR_ALIGN.out.log_final // channel: [ val(meta), path(log_final) ]
    star_align_log_out           = STAR_ALIGN.out.log_out // channel: [ val(meta), path(log_out) ]
    star_align_log_progress      = STAR_ALIGN.out.log_progress // channel: [ val(meta), path(log_progress) ]
    star_align_wig               = STAR_ALIGN.out.wig // channel: [ val(meta), path(wig) ]
    star_align_bedgraph          = STAR_ALIGN.out.bedgraph // channel: [ val(meta), path(bg) ]

    umitools_dedup_log           = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.deduplog // channel: [ val(meta), path(log) ]
    samtools_stats               = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.stats // channel: [ val(meta), path(stats) ]
    flagstat                     = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats                     = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    ch_bam_bai                   = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam.join(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bai) // channel: [ val(meta), path(bam), path(bai) ]
    ch_cram_crai                 = SAMTOOLS_CONVERT.out.cram.join(SAMTOOLS_CONVERT.out.crai) // channel: [ val(meta), path(cram), path(crai) ]

    versions                     = ch_versions // channel: [ versions.yml, versions.yml, ... ]
}
