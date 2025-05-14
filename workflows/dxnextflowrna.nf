/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULES, alphabetical order
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'

// SUBWORKFLOWS, alphabetical order
include { BAM_QUANTIFICATION_FEATURECOUNTS } from '../subworkflows/local/bam_quantification_featurecounts'
include { FASTQ_BAM_QC                  } from '../subworkflows/local/fastq_bam_qc'
include { FASTQ_TRIM_FILTER_ALIGN_DEDUP } from '../subworkflows/local/fastq_trim_filter_align_dedup'

// FUNCTIONS, alphabetical order
include { methodsDescriptionText        } from '../subworkflows/local/utils_umcugenetics_dxnextflowrna_pipeline'
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DXNEXTFLOWRNA {
    main:
    def createMetaWithIdSimpleName = { file -> [[id: file.getSimpleName()], file] }

    // Reference file channels
    ch_fasta_fai = Channel
        .fromPath(params.fasta)
        .combine(Channel.fromPath(params.fai))
        .map { fasta, fai -> [[id: fasta.getSimpleName()], fasta, fai] }
    ch_gene_bed = Channel.fromPath(params.gene_bed)
    ch_gtf = Channel
        .fromPath(params.gtf)
        .map(createMetaWithIdSimpleName)
        .first()
    ch_star_index = Channel
        .fromPath(params.star_index)
        .map(createMetaWithIdSimpleName)
        .first()
    ch_ref_flat = Channel.fromPath(params.ref_flat).first()
    ch_rrna_interval = Channel.fromPath(params.rrna_intervals).first()
    ch_sortmerna_fastas = Channel
        .fromPath(params.rrna_database_manifest)
        .splitCsv(by: 1, strip: true)
        .map { line -> file(line[0], checkIfExists: true) }
        .collect()
        .map { files -> [[id: file(params.rrna_database_manifest).getSimpleName()], files] }
    ch_sortmerna_index = Channel
        .fromPath(params.sortmerna_index)
        .map(createMetaWithIdSimpleName)
        .first()

    // Input channel
    ch_fastq = Channel
        .fromFilePairs("${params.input}/*_R{1,2}_001.fastq.gz")
        .map { meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            }
            else {
                fmeta.single_end = false
            }
            return [fmeta, fastq]
        }
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Run fastq_trim_filter_align_dedup
    //
    FASTQ_TRIM_FILTER_ALIGN_DEDUP(
        ch_fasta_fai,
        ch_fastq,
        ch_gtf,
        ch_sortmerna_fastas,
        ch_sortmerna_index,
        ch_star_index,
        params.seq_platform,
        params.seq_center,
        false,
    )
    ch_versions = ch_versions.mix(FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.versions)

    // Add fastq_trim_filter_align_dedup results to MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.trim_log.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.trim_zip.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.sortmerna_log.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.star_align_log_final.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.star_align_read_per_gene_tab.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.umitools_dedup_log.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.samtools_stats.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.flagstat.collect { it[1] }.ifEmpty([]),
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.idxstats.collect { it[1] }.ifEmpty([]),
    )

    //
    // SUBWORKFLOW: Run fastq_bam_qc
    //
    FASTQ_BAM_QC(
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.ch_bam_bai,
        ch_fasta_fai.map { meta, fasta, fai -> [fasta] },
        ch_fastq,
        ch_gene_bed,
        ch_ref_flat,
        ch_rrna_interval,
    )
    ch_versions = ch_versions.mix(FASTQ_BAM_QC.out.versions)

    // Add fastq_bam_qc results to MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQ_BAM_QC.out.fastqc_zip.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.bamstat_txt.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.inferexperiment_txt.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.innerdistance_freq.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.junctionannotation_all.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.junctionsaturation_rscript.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.readdistribution_txt.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.readduplication_pos_xls.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.tin_txt.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.rna_metrics.collect { it[1] }.ifEmpty([]),
        FASTQ_BAM_QC.out.lc_extrap.collect { it[1] }.ifEmpty([]),
    )

    //
    // SUBWORKFLOW: Run bam_quantification_featurecounts
    //
    BAM_QUANTIFICATION_FEATURECOUNTS(
        FASTQ_TRIM_FILTER_ALIGN_DEDUP.out.ch_bam_bai,
        ch_gtf
    )
    ch_versions = ch_versions.mix(BAM_QUANTIFICATION_FEATURECOUNTS.out.versions)





    // Add bam_quantification_featurecounts results to MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(
        BAM_QUANTIFICATION_FEATURECOUNTS.out.gene_counts_summary.collect { it[1] }.ifEmpty([]),
        BAM_QUANTIFICATION_FEATURECOUNTS.out.exon_counts_summary.collect { it[1] }.ifEmpty([]),
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    // Collect workflow params
    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    // Add methods description
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )
    // Collate software versions
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report // channel: [ /path/to/multiqc_report.html]
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
