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
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from './subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import local subworkflows, alphabetical order
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_BAM_QC } from './subworkflows/local/fastq_bam_qc'
include { FASTQ_TRIM_FILTER_ALIGN_DEDUP } from './subworkflows/local/fastq_trim_filter_align_dedup'

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
    ch_sortmerna_fastas = Channel
        .fromPath( params.rrna_database_manifest )
        .splitCsv( by: 1 , strip: true)
        .map { line -> file(line[0], checkIfExists: true) }
        .collect()
        .map { files -> [[id: file(params.rrna_database_manifest).getSimpleName()], files] }
    ch_sortmerna_index = Channel.fromPath(params.sortmerna_index)
        .map(createMetaWithIdSimpleName)
        .first()

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
    FASTQ_TRIM_FILTER_ALIGN_DEDUP(
        ch_fasta_fai,
        ch_fastq,
        ch_gtf,
        ch_sortmerna_fastas,
        ch_sortmerna_index,
        ch_star_index,
        params.seq_platform,
        params.seq_center,
        false
    )

    FASTQ_BAM_QC(
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.ch_bam_bai,
        ch_fasta_fai.map { meta, fasta, fai -> [ fasta ] },
        ch_fastq,
        ch_gene_bed,
        ch_ref_flat,
        ch_rrna_interval
    )

    // MultiQC
    // Collect workflow params
    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    // Collate software versions
    ch_collated_versions = softwareVersionsToYAML(
		Channel.empty().mix(
            Channel.fromPath(params.sortmerna_index_versions),
		    FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.versions,
	        FASTQ_BAM_QC.out.versions
		)
	)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    // Run MultiQC
    MULTIQC(
        Channel.empty().mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
	        ch_collated_versions,
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.trim_log.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.sortmerna_log.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.star_align_log_final.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.star_align_read_per_gene_tab.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.umitools_dedup_log.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.samtools_stats.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.flagstat.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.idxstats.collect{it[1]}.ifEmpty([]),
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
