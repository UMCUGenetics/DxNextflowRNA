#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq module
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// bam module
include extractBamFromDir from './NextflowModules/Utils/bam.nf'

// QC modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MultiQC as MultiQC_pre from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'-o preQC')
include MultiQC as MultiQC_post from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'-o postQC')

// Trimming module
include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params(optional: '--fastqc', single_end: false)

// rRNA removal module
include SortMeRNA from './NextflowModules/SortMeRNA/4.3.3/SortMeRNA.nf' params( single_end: false )

// Mapping modules
include GenomeGenerate from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf'
include AlignReads from './NextflowModules/STAR/2.7.3a/AlignReads.nf' params(optional: '--outReadsUnmapped Fastx', single_end: false)
include Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Flagstat as Flagstat_raw from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'
include BWAMapping from './NextflowModules/BWA-Mapping/bwa-0.7.17_samtools-1.9/Mapping.nf' params(
    genome_fasta: "$params.genome_fasta", optional: '-c 100 -M'
)
include Merge as Sambamba_Merge_RNA from './NextflowModules/Sambamba/0.7.0/Merge.nf'
include Merge as Sambamba_Merge_WES from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// After mapping QC
include RSeQC from './NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( single_end: false)
include LCExtrap from './NextflowModules/Preseq/2.0.3/LCExtrap.nf' params( optional:'-v -B -D')
include GtfToGenePred from './NextflowModules/UCSC/377/GtfToGenePred.nf'
include GenePredToBed from './NextflowModules/UCSC/377/GenePredToBed.nf'
include CollectMultipleMetrics as PICARD_CollectMultipleMetrics from './NextflowModules/Picard/2.22.0/CollectMultipleMetrics.nf' params(
    genome: "$params.genome_fasta",
    optional: "PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE"
)
include EstimateLibraryComplexity as PICARD_EstimateLibraryComplexity from './NextflowModules/Picard/2.22.0/EstimateLibraryComplexity.nf' params(
    optional: "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
)
include FixMateInformation as PICARD_FixMateInformation from './NextflowModules/Picard/2.22.0/FixMateInformation.nf' params(optional: '')
include UmiAwareMarkDuplicatesWithMateCigar as PICARD_UmiAwareMarkDuplicatesWithMateCigar from './NextflowModules/Picard/2.22.0/UmiAwareMarkDuplicatesWithMateCigar.nf' params(optional:'')
//include Qualimap from './NextflowModules/qualimap-2.2.2d.1.nf' 
include FastqScreen from './NextflowModules/Fastq-screen-0.15.3.nf'

// WES

// IndelRealignment modules
include RealignerTargetCreator as GATK_RealignerTargetCreator from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/RealignerTargetCreator.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta", optional: "$params.gatk_rtc_options"
)
include IndelRealigner as GATK_IndelRealigner from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/IndelRealigner.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta", optional: ""
)
include ViewUnmapped as Sambamba_ViewUnmapped from './NextflowModules/Sambamba/0.7.0/ViewUnmapped.nf'
include Merge as Sambamba_Merge from './NextflowModules/Sambamba/0.7.0/Merge.nf'

// HaplotypeCaller modules
include IntervalListTools as PICARD_IntervalListTools from './NextflowModules/Picard/2.22.0/IntervalListTools.nf' params(
    scatter_count: "500", optional: ""
)
include HaplotypeCaller as GATK_HaplotypeCaller from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/HaplotypeCaller.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta", optional: "$params.gatk_hc_options"
)
include VariantFiltrationSnpIndel as GATK_VariantFiltration from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/VariantFiltration.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta", snp_filter: "$params.gatk_snp_filter",
    snp_cluster: "$params.gatk_snp_cluster", indel_filter: "$params.gatk_indel_filter"
)
include CombineVariants as GATK_CombineVariants from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/CombineVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta", optional: "--assumeIdenticalSamples"
)
include SelectVariantsSample as GATK_SingleSampleVCF from './NextflowModules/GATK/3.8-1-0-gf15c1c3ef/SelectVariants.nf' params(
    gatk_path: "$params.gatk_path", genome: "$params.genome_fasta"
)
include Tabix from './NextflowModules/Tabix.nf'

// DROP
include DROP from './NextflowModules/drop-1.2.4.nf'

def analysis_id = params.outdir.split('/')[-1]
def fastq_files = extractFastqPairFromDir(params.fastq_path)
def wes_files = extractFastqPairFromDir(params.wes_path)

// Define chromosomes used to scatter GATK_RealignerTargetCreator
def chromosomes = Channel.fromPath(params.genome_fasta.replace('fasta', 'dict'))
    .splitCsv(sep:'\t', skip:1)
    .map{type, chr, chr_len, md5, file -> [chr.minus('SN:')]}

workflow {
    // Use bam files if bam is set true
    if (!params.bam) {
        // FASTQC
        FastQC(fastq_files)

        // MultiQC pre trimming
        pre_QC = Channel.empty().mix(FastQC.out).collect()
        MultiQC_pre(analysis_id, pre_QC)

        // TrimGalore
        TrimGalore(fastq_files)
        trimmed_fastqs = TrimGalore.out.fastqs_trimmed

        //Remove rRNA with SortMeRNA
        rRNA_database = file(params.rRNA_database_manifest)
        if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
        sortmerna_fasta = Channel
            .from( rRNA_database.readLines() )
            .map { row -> file(row) }
        SortMeRNA(trimmed_fastqs, sortmerna_fasta.collect())
        final_fastqs = SortMeRNA.out.non_rRNA_fastqs

        //Create STAR index if not present
        if (params.star_index) {
            star_index = Channel
                .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}"}
        } else {
            // Import genome
            genome_fasta = Channel
                .fromPath(params.genome_fasta, checkIfExists: true)
                .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}"}
            // Import GTF file
            genome_gtf = Channel
                .fromPath( params.genome_gtf, checkIfExists: true )
                .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}
        
            GenomeGenerate(genome_fasta, genome_gtf)
            star_index = GenomeGenerate.out.star_index
        }

        // Mapping, sorting, indexing and stats
        AlignReads( final_fastqs, star_index.collect(), genome_gtf.collect() )
        Index( AlignReads.out.bam_file.map {sample_id, rg_id, bam ->
                                           [sample_id, bam] })
        bam_files = AlignReads.out.bam_file.map {sample_id, rg_id, bam ->
                                           [sample_id, bam] }
        bam_files = bam_files.join(Index.out)
        
        //Sambamba_Merge_RNA(bam_files.groupTuple())
        //bam_files = Sambamba_Merge.out

        Flagstat_raw(bam_files)
        star_logs = AlignReads.out.log.mix(AlignReads.out.final_log)
        flagstat_logs = Flagstat_raw.out
    }
    else{
        bam_files = extractBamFromDir(params.fastq_path)
    }

    // Remove duplicates with Picard and stats
    PICARD_CollectMultipleMetrics(bam_files)
    PICARD_EstimateLibraryComplexity(bam_files)
    //PICARD_FixMateInformation(bam_files)
    //PICARD_UmiAwareMarkDuplicatesWithMateCigar(PICARD_FixMateInformation.out.fixmate_bam)

    // Load gtf file
    genome_gtf = Channel
        .fromPath( params.genome_gtf, checkIfExists: true )
        .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}

    // After mapping quality using RSeQC, Preseq (LCExtrap) and Qualimap
    if (params.genome_bed ) {
        //Create bed12 index file
        genome_bed = Channel
            .fromPath(params.genome_bed, checkIfExists: true)
            .ifEmpty { exit 1, "Genome bed12 file not found: ${params.genome_bed}"}
    } else if ( !params.genome_bed ) {
        GtfToGenePred(genome_gtf)
        GenePredToBed(GtfToGenePred.out.genome_genepred)
        genome_bed = GenePredToBed.out.genome_bed12
    }
    RSeQC(bam_files, genome_bed.collect())
    LCExtrap(bam_files)

    // Qualimap
    //Qualimap(bam_files, genome_gtf)

    // MultiQC
    
    if (!params.bam) {
        FastqScreen(final_fastqs, params.fastq_screen_config)
        qc_files = Channel.empty().mix(TrimGalore.out.trimming_report, TrimGalore.out.fastqc_report, SortMeRNA.out.qc_report, star_logs, flagstat_logs, RSeQC.out, LCExtrap.out, PICARD_CollectMultipleMetrics.out, PICARD_EstimateLibraryComplexity.out, FastqScreen.out).collect()
    }
    else{
        FastqScreen(fastq_files, params.fastq_screen_config)
        qc_files = Channel.empty().mix(PICARD_CollectMultipleMetrics.out, PICARD_EstimateLibraryComplexity.out, RSeQC.out, LCExtrap.out, FastqScreen.out).collect()
    }
    MultiQC_post(analysis_id, qc_files)

    // WES data
    BWAMapping(wes_files)
    wes_bam = BWAMapping.out.map{
            sample_id, rg_id, bam_file, bai_file -> [sample_id, bam_file, bai_file]}.groupTuple()

    GATK_RealignerTargetCreator(wes_bam.combine(chromosomes))
    GATK_IndelRealigner(wes_bam.combine(GATK_RealignerTargetCreator.out, by: 0))
    Sambamba_ViewUnmapped(wes_bam)
    Sambamba_Merge_WES(GATK_IndelRealigner.out.mix(Sambamba_ViewUnmapped.out).groupTuple())

    PICARD_IntervalListTools(Channel.fromPath("$params.gatk_hc_interval_list"))
    GATK_HaplotypeCaller(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [analysis_id, bam_file, bai_file]}
            .groupTuple()
            .combine(PICARD_IntervalListTools.out.flatten())
    )
    GATK_VariantFiltration(GATK_HaplotypeCaller.out)
    GATK_CombineVariants(GATK_VariantFiltration.out.groupTuple())
    GATK_SingleSampleVCF(GATK_CombineVariants.out.combine(
        Sambamba_Merge.out.map{sample_id, bam_file, bai_file -> [sample_id]})
    )
    Tabix(GATK_SingleSampleVCF.out)

    // DROP
    //DROP("/hpc/diaggen/projects/RNAseq_Jade/drop/Data", "/hpc/diaggen/projects/RNAseq_Jade/drop/Scripts", "/hpc/diaggen/projects/RNAseq_Jade/drop/.drop", "/hpc/diaggen/projects/RNAseq_Jade/drop/Snakefile", "/hpc/diaggen/projects/RNAseq_Jade/drop/config.yaml", "/hpc/diaggen/projects/RNAseq_Jade/drop/R_lib", "/hpc/diaggen/projects/RNAseq_Jade/drop/readme.md")

    // rMATS
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "RNAseq Workflow Successful: ${analysis_id}"
        sendMail(
            to: params.email.trim(),
            subject: subject,
            body: email_html,
            attach: ["${params.outdir}/QC/preQC/${analysis_id}_multiqc_report.html","${params.outdir}/QC/postQC/${analysis_id}_multiqc_report.html"]
        )

    } else {
        def subject = "RNAseq Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
