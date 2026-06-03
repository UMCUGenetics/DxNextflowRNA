include { GATK4_BEDTOINTERVALLIST        } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS        } from '../../../modules/nf-core/gatk4/intervallisttools/main'
include { BEDOPS_GTF2BED                 } from '../../../modules/nf-core/bedops/gtf2bed/main'

workflow PREPARE_REFERENCES {
    take:
    ch_fasta_fai
    ch_gtf

    main:

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })
    def ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    BEDOPS_GTF2BED(ch_gtf)
    GATK4_BEDTOINTERVALLIST(BEDOPS_GTF2BED.out.bed.collect(), ch_dict)

    // Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)



    emit:
    dict = ch_dict
    bed = BEDOPS_GTF2BED.out.bed.collect()
    interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map { _meta, bed -> [bed] }.collect()

}
