/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULES, alphabetical order
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE } from '../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_EXON } from '../../modules/nf-core/subread/featurecounts/main'

workflow BAM_QUANTIFICATION_FEATURECOUNTS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ] ]
    ch_gtf // channel: [ val(meta), path(gtf) ]

    main:
    // Create empty versions channel, and fill with each tools version
    ch_versions = Channel.empty()
    path_gtf = ch_gtf.map { meta, gtf -> [gtf] }
    ch_bam_gtf = ch_bam_bai.map { meta, bam, bai -> [meta, bam, path_gtf] }

    // Run subread featurecounts with different meta_features.
    SUBREAD_FEATURECOUNTS_GENE(ch_bam_gtf)
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_GENE.out.versions.first())

    SUBREAD_FEATURECOUNTS_EXON(ch_bam_gtf)
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_EXON.out.versions.first())

    emit:
    gene_counts_summary = SUBREAD_FEATURECOUNTS_GENE.out.summary // path featureCounts.txt.summary
    gene_counts         = SUBREAD_FEATURECOUNTS_GENE.out.counts // channel: [ val(meta), path(featureCounts.tsv) ]
    exon_counts_summary = SUBREAD_FEATURECOUNTS_EXON.out.summary // path featureCounts.txt.summary
    exon_counts         = SUBREAD_FEATURECOUNTS_EXON.out.counts // channel: [ val(meta), path(featureCounts.tsv) ]
    versions            = ch_versions // channel: [ versions.yml ]
}
