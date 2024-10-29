#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowRNA
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import local subworkflows, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_BAM_QC } from './subworkflows/local/fastq_bam_qc'
include { FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA } from './subworkflows/local/FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    def createMetaWithIdSimpleName = {file -> [[id: file.getSimpleName()], file]}

    // Reference file channels
    ch_fasta_fai = Channel.fromPath(params.fasta)
        .combine(Channel.fromPath(params.fai))
        .map{fasta, fai -> [[id: fasta.getSimpleName()], fasta, fai]}
    ch_gene_bed = Channel.fromPath(params.gene_bed)
    ch_gtf = Channel.fromPath(params.gtf)
        .map(createMetaWithIdSimpleName)
        .first()
    ch_star_index = Channel.fromPath(params.star_index)
        .map(createMetaWithIdSimpleName)
        .first()
    ch_ref_flat = Channel.fromPath(params.ref_flat).first()
    ch_rrna_interval = Channel.fromPath(params.rrna_intervals).first()

    // Input channel
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
            return [ fmeta, fastq ]
        }

    // Subworkflows
    FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA(
        ch_fasta_fai, ch_fastq, ch_gtf, ch_star_index, params.seq_platform, params.seq_center, false
    )
    FASTQ_BAM_QC(
        FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.ch_bam_bai,
        ch_fasta_fai.map { meta, fasta, fai -> [ fasta ] },
        ch_fastq,
        ch_gene_bed,
        ch_ref_flat,
        ch_rrna_interval
    )

    // MultiQC
    MULTIQC(
        Channel.empty().mix(
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.versions,
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.trim_log.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.sortmerna_log.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.star_align_log_final.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_TRIMGALORE_STAR_SORTMERNA.out.star_align_read_per_gene_tab.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.versions,
            FASTQ_BAM_QC.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.bamstat_txt.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.inferexperiment_txt.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.innerdistance_freq.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.junctionannotation_all.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.junctionsaturation_rscript.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.readdistribution_txt.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.readduplication_pos_xls.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.tin_txt.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.rna_metrics.collect{it[1]}.ifEmpty([]),
            FASTQ_BAM_QC.out.lc_extrap.collect{it[1]}.ifEmpty([]),
        ).collect(),
        Channel.fromPath("${params.multiqc_yaml}", checkIfExists: true),
        [],  // extra_multiqc_config
        [],  // multiqc_logo
        [],  // replace_names
        []   // sample_names
    )
}
