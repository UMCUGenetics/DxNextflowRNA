/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FASTQ_BAM_QC {
    take:
    ch_bam  // channel: [ val(meta), [ bam ] ]
    ch_fasta  // channel: path(fasta)
    ch_fastq  // channel: [ val(meta), [ fastq_forward, fastq_reverse ] ]
    ch_gene_bed  // channel: path(gene_bed.bed)
    ch_ref_flat  // channel: path(refFlat.txt)
    ch_rrna_interval  // channel: path(rrna_intervals.bed)

    main:
    ch_versions = Channel.empty()

    FASTQC(ch_fastq)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    PICARD_COLLECTRNASEQMETRICS(ch_bam, ch_ref_flat, ch_fasta, ch_rrna_interval)
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    fastqc_html = FASTQC.out.html  // channel: [ val(meta), [ html ] ]
    fastqc_zip = FASTQC.out.zip  // channel: [ val(meta), [ zip ] ]

    rna_metrics = PICARD_COLLECTRNASEQMETRICS.out.metrics  // channel: [ val(meta), [ rna_metrics ]]

    versions = ch_versions  // channel: [ versions.yml ]
}
