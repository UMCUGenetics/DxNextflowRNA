#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowRNA
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from './modules/nf-core/fastqc/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'

include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS } from './subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'

include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'

include { SUBREAD_FEATURECOUNTS } from './modules/nf-core/subread/featurecounts/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
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
            [ fmeta, fastq ]
        }

    //TRIMGALORE

    FASTQC(ch_fastq)

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )

    ch_star_index = file('/hpc/diaggen/data/databases/STAR_ref_genome_index/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_2.7.9a')
    ch_star_index = ch_star_index.map {
        star_index -> [star_index.getSimpleName(), star_index]
    }

    ch_gtf = file('/hpc/diaggen/data/databases/gencode/gencode.v44.primary_assembly.basic.annotation.gtf')
    ch_gtf = ch_gtf.map {
        gtf -> [gtf.getSimpleName(), gtf]
    }

    STAR_ALIGN (
        ch_fastq,
        ch_star_index,
        ch_gtf,
        false,
        'illumina',
        'UMCU Genetics'
    )

    SAMTOOLS_MERGE( STAR_ALIGN.out.bam_sorted )

    SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set { ch_bam_bai }

    BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS(
        ch_bam_bai,
        true
    )

    SUBREAD_FEATURECOUNTS(
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS.out.bam.map{
        meta, bam -> [ meta, bam, params.genome ]
        }
    )

}
