#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq modules
//include extractFastqPairFromDir from './modules/fastq.nf'
//def fastq_files = extractFastqPairFromDir(params.fastq_path)

//QC module
//include FastQC from './modules/FastQC.nf' params(optional:'')

//bcftools view multiallelic  ## ME only used for monoallelic expresion
include norm_bcf_norm as BCFtools_Norm_m_both from './modules/bcftools.nf' params(optional: "-m-both")
include norm_bcf_norm as BCFtools_Norm_d_both from './modules/bcftools.nf' params(optional: "-d both")
include view_bcf_norm as BCFtools_view from './modules/bcftools_view.nf' params(optional: "-m2 -M2 -v snps")

//samtools bgzip ## ME only used for monoallelic expresion
include bgzip from './modules/samtools_bgzip.nf' params(optional: "-c")
include tabix from './modules/samtools_tabix.nf'

//GATK ## ME only used for monoallelic expresion
include ASEReadCounter from './modules/ASEReadCounter.nf'

//MAE ## ME only used for monoallelic expresion
include MAE from './modules/MAE.nf'

//ASGL ## ME for variantcalling
include ASGAL from './modules/ASGAL.nf' params(optional: "--multi")

transcriptome_ch = Channel.fromPath(params.transcripts)
genome_ch = Channel.fromPath(params.genome)
gtf_ch = Channel.fromPath(params.gtf)
//STAR_genome_index_ch = Channel.fromPath(params.STAR_genome_index)
//vcf_file_ch=Channel.fromPath(params.vcf_infile)
//bam_infile_ch=Channel.fromPath(params.bam_infile)

workflow {

    // Monoallelic expresion
    norm_bcf_norm_m_ch=BCFtools_Norm_m_both(params.vcf_infile)
    norm_bcf_norm_d_ch=BCFtools_Norm_d_both(norm_bcf_norm_m_ch)
    bcf_view_ch=BCFtools_view(norm_bcf_norm_d_ch)
    bgzip_ch=bgzip(bcf_view_ch)
    tabix_ch=tabix(bgzip_ch)
    count_ch=ASEReadCounter(genome_ch, params.vcf_infile, params.bam_infile)
    MAE_ch=MAE(count_ch)

    //ASGAL_ch=ASGAL(genome_ch, transcriptome_ch, gtf_ch, fastq_files)
}
