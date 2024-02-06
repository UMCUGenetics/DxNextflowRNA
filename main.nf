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
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { OUTRIDER as OUTRIDER_EXON } from './NextflowModules/Outrider/1.20.0/main'
include { OUTRIDER as OUTRIDER_GENE } from './NextflowModules/Outrider/1.20.0/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_EXON } from './modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE } from './modules/nf-core/subread/featurecounts/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Reference file channels
    ch_star_index = Channel.fromPath(params.star_index).map {star_index -> [star_index.getSimpleName(), star_index] }
    ch_gtf = Channel.fromPath(params.gtf).map { gtf -> [gtf.getSimpleName(), gtf] }

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
            [ fmeta, fastq ]
        }

    // Trimgalore
    TRIMGALORE(
        ch_fastq
    )
    
    //Star Align
    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index.first(),
        ch_gtf.first(),
        false,
        'illumina',
        'UMCU Genetics'
    )

    //Samtools
    SAMTOOLS_MERGE(
        STAR_ALIGN.out.bam_sorted.map {
            meta, bam ->
                new_id = meta.id.split('_')[0]
                [ meta + [id: new_id], bam ]
        }.groupTuple(),
        [ [ id:'null' ], []],
        [ [ id:'null' ], []],
    )

    SAMTOOLS_INDEX ( SAMTOOLS_MERGE.out.bam )

    SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set { ch_bam_bai }

    //Featurecounts
     SUBREAD_FEATURECOUNTS_GENE(
        SAMTOOLS_MERGE.out.bam.map{
            meta, bam -> [ meta, bam, params.gtf, 'gene_id' ]
        }
    )

    SUBREAD_FEATURECOUNTS_EXON(
        SAMTOOLS_MERGE.out.bam.map{
            meta, bam -> [ meta, bam, params.gtf, 'exon_id' ]
        }
    )

    ch_outrider_ref_gene = params.refgene.contains(",") ? Channel.fromPath(params.refgene?.split(',') as List).view() : Channel.fromPath("$params.refgene/*.txt").view()
    ch_outrider_ref_exon = params.refexon.contains(",") ? Channel.fromPath(params.refexon?.split(',') as List).view() : Channel.fromPath("$params.refexon/*.txt").view()

    // Outrider Gene
    OUTRIDER_GENE(
        SUBREAD_FEATURECOUNTS_GENE.out.counts,
        ch_outrider_ref_gene.collect(),
        'gene'
    )

    // Outrider Gene
    OUTRIDER_EXON(
        SUBREAD_FEATURECOUNTS_EXON.out.counts,
        ch_outrider_ref_exon.toList(),
        'exon'
    )

    // QC
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

}