#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Indexing modules
include { INDEX } from './modules/salmon_index.nf'

// Quantification modules
include { QUANT } from './modules/salmon_quant.nf'

// fastq modules
include extractFastqPairFromDir from './modules/fastq.nf'

def fastq_files = extractFastqPairFromDir(params.fastq_path)

//QC module
include FastQC from './modules/FastQC.nf' params(optional:'')

//build star index
include GenomeGenerate from './modules/GenomeGenerate.nf'

//mapping using star
include AlignReads from './modules/AlignReads.nf'

//AS calling_rmats
include RMATS from './modules/rmats.nf'

//bcftools view multiallelic
include norm_bcf_norm as BCFtools_Norm_m_both from './modules/bcftools.nf' params(optional: "-m-both")
include norm_bcf_norm as BCFtools_Norm_d_both from './modules/bcftools.nf' params(optional: "-d both")
include view_bcf_norm as BCFtools_view from './modules/bcftools_view.nf' params(optional: "-m2 -M2 -v snps")

//samtools bgzip
include bgzip from './modules/samtools_bgzip.nf' params(optional: "-c")
include tabix from './modules/samtools_tabix.nf'


//GATK
include ASEReadCounter from './modules/ASEReadCounter.nf'

//MAE
include MAE from './modules/MAE.nf'

//ASGL
include ASGAL from './modules/ASGAL.nf' params(optional: "--multi")


transcriptome_ch = Channel.fromPath(params.transcripts)
genome_ch = Channel.fromPath(params.genome)
gtf_ch = Channel.fromPath(params.gtf)
s1_ch = Channel.fromPath(params.s1)
s2_ch = Channel.fromPath(params.s2)
STAR_genome_index_ch = Channel.fromPath(params.STAR_genome_index)
vcf_file_ch=Channel.fromPath(params.vcf_infile)
bam_infile_ch=Channel.fromPath(params.bam_infile)

workflow {


// run
    fastqc_ch=FastQC(fastq_files)
    star_index_ch=GenomeGenerate(genome_ch,gtf_ch)
    star_mapping_ch=AlignReads(fastq_files,star_index_ch,gtf_ch)
    rmats_AS_calling_ch=RMATS(s1_ch, s2_ch, gtf_ch, STAR_genome_index_ch)
    index_ch=INDEX(transcriptome_ch)
    quant_ch=QUANT(index_ch,fastq_files)
    norm_bcf_norm_m_ch=BCFtools_Norm_m_both(vcf_file_ch)
    norm_bcf_norm_d_ch=BCFtools_Norm_d_both(norm_bcf_norm_m_ch)
    bcf_view_ch=BCFtools_view(norm_bcf_norm_d_ch)
    bgzip_ch=bgzip(bcf_view_ch)
    tabix_ch=tabix(bgzip_ch)
    count_ch=ASEReadCounter(genome_ch, vcf_file_ch, bam_infile_ch)
    MAE_ch=MAE(count_ch)
    ASGAL_ch=ASGAL(genome_ch, transcriptome_ch, gtf_ch, fastq_files)
}
