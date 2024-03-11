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
    def create_meta_with_id_name = {file -> [[id: file.getSimpleName()], file]}

    // Reference path channels as value channels
    ch_star_index = Channel.fromPath(params.star_index)
        .map(create_meta_with_id_name)
        .first()

    ch_gtf = Channel.fromPath(params.gtf)
        .map(create_meta_with_id_name)
        .first()
    
    ch_bed = Channel.fromPath(params.bed_file)
        .map(create_meta_with_id_name)
        .first()

    ch_genome_fasta = Channel.fromPath(params.genome)
        .map(create_meta_with_id_name)
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
            [ fmeta, fastq ]
        }
    
    // apply trimming to fastq files  
    TRIMGALORE(ch_fastq)
    
    // align reads to reference genome
    STAR_ALIGN(
        TRIMGALORE.out.reads,
        ch_star_index,
        ch_gtf,
        false,
        'illumina',
        'UMCU Genetics'
    )

    // merge all lanes for a sample
    // TODO: replace input with fasta and fai instead of empty.
    SAMTOOLS_MERGE(
        STAR_ALIGN.out.bam.map {
            meta, bam ->
                new_id = meta.id.split('_')[0]
                [ meta + [id: new_id], bam ]
        }.groupTuple(),
        [ [ id:'null' ], []],
        [ [ id:'null' ], []],
    )

    // generate bai file
    BAM_SORT_STATS_SAMTOOLS(SAMTOOLS_MERGE.out.bam, ch_genome_fasta)

    // ch_bam_bai channel: [ val(meta), [ bam, bai ] ]
    BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.bai)
        .set { ch_bam_bai }
    

    // PreSeq LCExtrap generate library complexity plot
    PRESEQ_LCEXTRAP(BAM_SORT_STATS_SAMTOOLS.out.bam)

    // QC
    FASTQC(ch_fastq)

    // bed file with with gene_bed created on GTF.
    BAM_RSEQC(
        ch_bam_bai,
        ch_bed,
        ['bam_stat', 'inferexperiment', 'innerdistance', 'junctionannotation', 'junctionsaturation', 'readdistribution', 'readduplication'] // 'tin'
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
