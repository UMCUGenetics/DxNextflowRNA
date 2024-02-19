#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowRNA
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { OUTRIDER } from './NextflowModules/Outrider/1.20.0/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'
include { SUBREAD_FEATURECOUNTS } from './modules/nf-core/subread/featurecounts/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    fastq_to_bam(ch_fastq)
    featurecounts(fastq_to_bam.out)
    outrider(featurecounts.out)
    QC(ch_fastq)
}


workflow fastq_to_bam{
    take: ch_fastq
    main:
    // Reference file channels
    ch_star_index = Channel.fromPath(params.star_index).map {star_index -> [star_index.getSimpleName(), star_index] }
    ch_gtf = Channel.fromPath(params.gtf).map { gtf -> [gtf.getSimpleName(), gtf] }
    
    TRIMGALORE(
        ch_fastq
    )
    
    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index.first(),
        ch_gtf.first(),
        false,
        'illumina',
        'UMCU Genetics'
    )

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

    emit:
	SAMTOOLS_MERGE.out.bam
}


workflow featurecounts{
    take: 
	samtools_merge_bam
    main: 
    samtools_merge_bam_gene = samtools_merge_bam.map { meta, bam -> [meta, bam, params.gtf, "gene"] }
    samtools_merge_bam_exon = samtools_merge_bam.map { meta, bam -> [meta, bam, params.gtf, "exon"] }

    SUBREAD_FEATURECOUNTS(
	samtools_merge_bam_gene.concat(samtools_merge_bam_exon)
    )

    emit: 
	SUBREAD_FEATURECOUNTS.out.counts
}


workflow outrider{
    take: 
        featurecounts

    main:
    ch_outrider_ref_gene = params.refgene.contains(",") ? Channel.fromPath(params.refgene?.split(',') as List).view() : Channel.fromPath("$params.refgene").view()
    ch_outrider_ref_exon = params.refexon.contains(",") ? Channel.fromPath(params.refexon?.split(',') as List).view() : Channel.fromPath("$params.refexon").view()

    OUTRIDER(
        featurecounts,
        ch_outrider_ref_gene.collect(),
        ch_outrider_ref_exon.collect()
    )
}


workflow featurecounts_entry{
    samtools_merge_bam = Channel.fromPath("$params.outdir/bam_files/*.bam").map { bam -> [[id:bam.getSimpleName()], bam] }

    featurecounts(samtools_merge_bam)
}


workflow outrider_entry{
    featurecounts_gene = Channel.fromPath("$params.outdir/feature_counts/*gene*.txt").map { counts -> [[id:counts.getSimpleName()], counts, "gene"] }
    featurecounts_exon = Channel.fromPath("$params.outdir/feature_counts/*exon*.txt").map { counts -> [[id:counts.getSimpleName()], counts, "exon"] }

    outrider(featurecounts_gene.concat(featurecounts_exon))
}


workflow QC{
    take: ch_fastq
    main:
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