include { GATK4_BEDTOINTERVALLIST        } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS        } from '../../../modules/nf-core/gatk4/intervallisttools/main'
include { BEDOPS_GTF2BED                 } from '../../../modules/nf-core/bedops/gtf2bed/main'

workflow PREPARE_REFERENCES {
    main:
    def createMetaWithIdSimpleName = { file -> [[id: file.getSimpleName()], file] }

    ch_fasta_fai = Channel
        .fromPath(params.fasta)
        .combine(Channel.fromPath(params.fai))
        .map { fasta, fai -> [[id: fasta.getSimpleName()], fasta, fai] }
        .collect()

    ch_gene_bed = Channel
        .fromPath(params.gene_bed)
        .collect()
    ch_gtf = Channel
        .fromPath(params.gtf)
        .map(createMetaWithIdSimpleName)
        .first()
        .collect()
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

    ch_dbsnp = channel
        .fromPath(params.dbsnp)
        .map(createMetaWithIdSimpleName)
        .first()
    ch_dbsnp_tbi = channel
        .fromPath("${params.dbsnp}.tbi")
        .map(createMetaWithIdSimpleName)
        .first()


    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })
    def ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    BEDOPS_GTF2BED(ch_gtf)
    GATK4_BEDTOINTERVALLIST(BEDOPS_GTF2BED.out.bed.collect(), ch_dict)

    // Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)


    emit:
    ch_fasta_fai = ch_fasta_fai
    ch_gene_bed = ch_gene_bed
    ch_gtf = ch_gtf
    ch_star_index = ch_star_index
    ch_ref_flat = ch_ref_flat
    ch_rrna_interval = ch_rrna_interval
    ch_sortmerna_fastas = ch_sortmerna_fastas
    ch_sortmerna_index = ch_sortmerna_index
    ch_dbsnp = ch_dbsnp
    ch_dbsnp_tbi = ch_dbsnp_tbi
    dict = ch_dict
    bed = BEDOPS_GTF2BED.out.bed.collect()
    interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map { _meta, bed -> [bed] }.collect()

}
