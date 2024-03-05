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
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'
include { SUBREAD_FEATURECOUNTS } from './modules/nf-core/subread/featurecounts/main'

include { PRESEQ_LCEXTRAP } from './modules/nf-core/preseq/lcextrap/main'
include { BAM_RSEQC } from './subworkflows/nf-core/bam_rseqc/main'
include { BAM_SORT_STATS_SAMTOOLS } from './subworkflows/nf-core/bam_sort_stats_samtools/main'                                                                                                                               

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Reference path channels as value channels
    ch_star_index = Channel.fromPath(params.star_index)
        .map { star_index -> [star_index.getSimpleName(), star_index] }
        .first()

    ch_gtf = Channel.fromPath(params.gtf)
        .map{ gtf -> [gtf.getSimpleName(), gtf] }
        .first

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
    
    // apply trimming to fastq files  
    TRIMGALORE(ch_fastq)
    
    // align reads to reference genome
    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index.first(),
        ch_gtf.first(),
        false,
        'illumina',
        'UMCU Genetics'
    )

    // merge all lanes for a sample
    // TODO: replace input with fasta and fai instead of empty.
    SAMTOOLS_MERGE(
        STAR_ALIGN.out.bam_sorted.map {
            meta, bam ->
                new_id = meta.id.split('_')[0]
                [ meta + [id: new_id], bam ]
        }.groupTuple(),
        [ [ id:'null' ], []],
        [ [ id:'null' ], []],
    )

    // generate bai file
    BAM_SORT_STATS_SAMTOOLS(SAMTOOLS_MERGE.out.bam, ch_genome_fasta)

    // joint bam and bai file
    SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set { ch_bam_bai }

    // sort bam file by position
    SAMTOOLS_SORT(
        // STAR_ALIGN.out.bam
        ch_bam_bai
    )

    // get statistics on bam file for preseq
    SAMTOOLS_STATS(
        ch_bam_bai
    )

    // get statistics on sorted bam file for preseq
    SAMTOOLS_STATS(
        SAMTOOLS_SORT.out.bam
    )


    // PreSeq LCExtrap generate library complexity plot
    PRESEQ_LCEXTRAP(
        // STAR_ALIGN.out.bam_sorted
        // SAMTOOLS_MERGE.out.bam
        SAMTOOLS_SORT.out.bam
    )

    // QC
    FASTQC(ch_fastq)

    // TODO: replace bed file with gene_bed created on GTF.
    BAM_RSEQC(
        ch_bam_bai,
	    Channel.fromPath("/hpc/diaggen/users/ellen/rnaseq_rseqc/2023-803CHX.bed").first(),
        ['bam_stat']
    )

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
        TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]),
        PRESEQ_LCEXTRAP.out.lc_extrap.collect{it[1]}.ifEmpty([])
    ).collect()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )

}
