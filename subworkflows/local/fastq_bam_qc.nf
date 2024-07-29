/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from '../../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FASTQ_BAM_QC {
    take:
    ch_fastq // channel: [ val(meta), [ path()fastq1), path(fastq1) ] ]

    main:
    ch_versions = Channel.empty()

    FASTQC(ch_fastq)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    fastqc_html = FASTQC.out.html  // channel: [ val(meta), [ html ] ]
    fastqc_zip = FASTQC.out.zip  // channel: [ val(meta), [ zip ] ]

    versions = ch_versions  // channel: [ versions.yml ]
}
