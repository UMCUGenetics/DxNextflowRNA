include { BCFTOOLS_INDEX                 } from '../../../modules/nf-core/bcftools/index/main'
include { GATK4_HAPLOTYPECALLER          } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_SPLITNCIGARREADS         } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_VARIANTFILTRATION        } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { SAMTOOLS_MERGE                 } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index/main'


// workflow adapted from github.com/nf-core/rnavar
workflow BAM_VARIANT_CALLING {
    take:
    ch_bam_bai
    ch_fasta_fai
    ch_dict
    interval_list_split
    ch_dbsnp
    ch_dbsnp_tbi


    main:
    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] }
    ch_fai   = ch_fasta_fai.map { meta, _fasta, fai -> [meta, fai] }


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

    def bam_splitncigar_interval = GATK4_SPLITNCIGARREADS.out.bam
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

    def splitncigar_merged_bam = SAMTOOLS_MERGE.out.bam

    SAMTOOLS_INDEX(splitncigar_merged_bam)

    def splitncigar_merged_bam_indices = SAMTOOLS_INDEX.out.index

    def splitncigar_merged_bam_bai = splitncigar_merged_bam.join(splitncigar_merged_bam_indices, failOnDuplicate: true, failOnMismatch: true)

    def haplotypecaller_interval_bam = splitncigar_merged_bam_bai
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

    GATK4_MERGEVCFS(
        ch_vcf_tbi.map { meta, vcfs, _tbis -> [meta, vcfs] },
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
