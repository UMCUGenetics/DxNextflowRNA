/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import subworkflows, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAM_RSEQC } from '../../subworkflows/nf-core/bam_rseqc/main'

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

    BAM_RSEQC(ch_bam_bai, ch_gene_bed, params.rseqc_modules)
    ch_versions = ch_versions.mix(BAM_RSEQC.out.versions.first())


    emit:
    // FASTQC
    // TODO nf-core: edit emitted channels
    fastqc_html = FASTQC.out.html  // channel: [ val(meta), [ html ] ]
    fastqc_zip = FASTQC.out.zip  // channel: [ val(meta), [ zip ] ]

    // PICARD_COLLECTRNASEQMETRICS
    rna_metrics = PICARD_COLLECTRNASEQMETRICS.out.metrics  // channel: [ val(meta), [ rna_metrics ]]

    // BAM_RSEQC
    bamstat_txt = BAM_RSEQC.out.bamstat_txt  // channel: [ val(meta), txt ]

    innerdistance_all = BAM_RSEQC.out.innerdistance_all  // channel: [ val(meta), {txt, pdf, r} ]
    innerdistance_distance = BAM_RSEQC.out.innerdistance_distance  // channel: [ val(meta), txt ]
    innerdistance_freq = BAM_RSEQC.out.innerdistance_freq  // channel: [ val(meta), txt ]
    innerdistance_mean = BAM_RSEQC.out.innerdistance_mean  // channel: [ val(meta), txt ]
    innerdistance_pdf = BAM_RSEQC.out.innerdistance_pdf  // channel: [ val(meta), pdf ]
    innerdistance_rscript = BAM_RSEQC.out.innerdistance_rscript  // channel: [ val(meta), r   ]

    inferexperiment_txt = BAM_RSEQC.out.inferexperiment_txt  // channel: [ val(meta), txt ]

    junctionannotation_all = BAM_RSEQC.out.junctionannotation_all  // channel: [ val(meta), {bed, xls, pdf, r, log} ]
    junctionannotation_bed = BAM_RSEQC.out.junctionannotation_bed  // channel: [ val(meta), bed ]
    junctionannotation_interact_bed = BAM_RSEQC.out.junctionannotation_interact_bed  // channel: [ val(meta), bed ]
    junctionannotation_xls = BAM_RSEQC.out.junctionannotation_xls  // channel: [ val(meta), xls ]
    junctionannotation_pdf = BAM_RSEQC.out.junctionannotation_pdf  // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf = BAM_RSEQC.out.junctionannotation_events_pdf  // channel: [ val(meta), pdf ]
    junctionannotation_rscript = BAM_RSEQC.out.junctionannotation_rscript  // channel: [ val(meta), r   ]
    junctionannotation_log = BAM_RSEQC.out.junctionannotation_log  // channel: [ val(meta), log ]

    junctionsaturation_all = BAM_RSEQC.out.junctionsaturation_all  // channel: [ val(meta), {pdf, r} ]
    junctionsaturation_pdf = BAM_RSEQC.out.junctionsaturation_pdf  // channel: [ val(meta), pdf ]
    junctionsaturation_rscript = BAM_RSEQC.out.junctionsaturation_rscript  // channel: [ val(meta), r   ]

    readdistribution_txt = BAM_RSEQC.out.readdistribution_txt  // channel: [ val(meta), txt ]

    readduplication_all = BAM_RSEQC.out.readduplication_all  // channel: [ val(meta), {xls, pdf, r} ]
    readduplication_seq_xls = BAM_RSEQC.out.readduplication_seq_xls  // channel: [ val(meta), xls ]
    readduplication_pos_xls = BAM_RSEQC.out.readduplication_pos_xls  // channel: [ val(meta), xls ]
    readduplication_pdf = BAM_RSEQC.out.readduplication_pdf  // channel: [ val(meta), pdf ]
    readduplication_rscript = BAM_RSEQC.out.readduplication_rscript  // channel: [ val(meta), r   ]

    tin_txt = BAM_RSEQC.out.tin_txt  // channel: [ val(meta), txt ]

    // Subworkflow specific outputs
    versions = ch_versions  // channel: [ versions.yml ]
}
