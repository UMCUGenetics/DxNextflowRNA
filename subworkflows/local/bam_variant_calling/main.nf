include { BCFTOOLS_INDEX                 } from '../../../modules/nf-core/bcftools/index/main'
include { GATK4_BEDTOINTERVALLIST        } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_INTERVALLISTTOOLS        } from '../../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_MERGEVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_SPLITNCIGARREADS         } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_VARIANTFILTRATION        } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GTF2BED                        } from '../../../modules/local/gtf2bed/main'
include { SAMTOOLS_MERGE                 } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'


// workflow adapted from github.com/nf-core/rnavar
workflow BAM_VARIANT_CALLING {
    take:
    ch_bam_bai
    ch_fasta_fai
    ch_gtf
    ch_dbsnp
    ch_dbsnp_tbi


    main:
    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] }
    ch_fai   = ch_fasta_fai.map { meta, _fasta, fai -> [meta, fai] }


    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })
    def ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    GTF2BED(ch_gtf, 'exon')

    GATK4_BEDTOINTERVALLIST(GTF2BED.out.bed.collect(), ch_dict)
    def interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list

    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
    GATK4_INTERVALLISTTOOLS(interval_list)

    def interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map { _meta, bed -> [bed] }.collect()


    def bam_interval = ch_bam_bai
        .combine(interval_list_split)
        .map { meta, bam, bai, intervals ->
            def new_meta = meta + [interval_count: intervals.size()]
            [new_meta, bam, bai, intervals]
        }
        .transpose(by: 3)
        .map { meta, bam, bai, interval ->
            [meta + [id: "${meta.id}_${interval.baseName}", sample: meta.id], bam, bai, interval]
        }

    GATK4_SPLITNCIGARREADS(
        bam_interval,
        ch_fasta,
        ch_fai,
        ch_dict,
    )

    def bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam

    def bam_splitncigar_interval = bam_splitncigar
        .map { meta, bam ->
            def new_meta = meta + [id: meta.sample] - meta.subMap('sample') - meta.subMap('interval_count')
            [groupKey(new_meta, meta.interval_count), bam]
        }
        .groupTuple()

    SAMTOOLS_MERGE(
        bam_splitncigar_interval,
        ch_fasta,
        ch_fai
    )

    def splitncigar_bam = SAMTOOLS_MERGE.out.bam

    SAMTOOLS_INDEX(splitncigar_bam)

    def splitncigar_bam_indices = SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.csi)
        .mix(SAMTOOLS_INDEX.out.crai)

    def splitncigar_bam_bai = splitncigar_bam.join(splitncigar_bam_indices, failOnDuplicate: true, failOnMismatch: true)

    def haplotypecaller_interval_bam = splitncigar_bam_bai
        .combine(interval_list_split)
        .map { meta, bam, bai, interval_lists ->
            [meta + [interval_count: interval_lists.size()], bam, bai, interval_lists]
        }
        .transpose(by: 3)
        .map { meta, bam, bai, interval_list_ ->
            [meta + [id: meta.id + "_" + interval_list_.baseName, sample: meta.id, variantcaller: 'haplotypecaller'], bam, bai, interval_list_, []]
        }

        GATK4_HAPLOTYPECALLER(
            haplotypecaller_interval_bam,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi,
        )

        def ch_vcf_tbi = GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true, failOnDuplicate: true)
            .map { meta, vcf, tbi ->
                [groupKey(meta + [id: meta.sample] - meta.subMap('sample', "interval_count"), meta.interval_count), vcf, tbi]
            }
        .groupTuple()


        def haplotypecaller_raw = ch_vcf_tbi.map { meta, vcfs, _tbis -> [meta, vcfs] }
        GATK4_MERGEVCFS(
            haplotypecaller_raw,
            ch_dict,
        )

    BCFTOOLS_INDEX(GATK4_MERGEVCFS.out.vcf)

    GATK4_VARIANTFILTRATION(
        GATK4_MERGEVCFS.out.vcf.join(BCFTOOLS_INDEX.out.tbi),
        ch_fasta,
        ch_fai,
        ch_dict,
        [[],[]]
    )




    emit:
    ch_vcf_tbi_unfiltered = ch_vcf_tbi
    ch_vcf_filtered = GATK4_VARIANTFILTRATION.out.vcf
    ch_tbi_filtered = GATK4_VARIANTFILTRATION.out.tbi

}
