#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq module
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// bam module
include extractBamFromDir from './NextflowModules/Utils/bam.nf'

// QC modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MultiQC as MultiQC_pre from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'')
include MultiQC as MultiQC_post from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'')

// Trimming module
include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params(optional:'--fastqc', single_end: params.single_end)

// rRNA removal module
include SortMeRNA from './NextflowModules/SortMeRNA/4.3.3/SortMeRNA.nf' params( single_end: params.single_end )

// Mapping modules
include GenomeGenerate from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf'
include AlignReads from './NextflowModules/STAR/2.7.3a/AlignReads.nf'
include Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Flagstat as Flagstat_raw from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'

// After mapping QC
include RSeQC from './NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( single_end:params.single_end)
include RSeQC_TIN from './NextflowModules/RSeQC/3.0.1/RSeQC.nf'
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
include Qualimap from './NextflowModules/qualimap-2.2.2d.1.nf' 
include FastqScreen from './NextflowModules/Fastq-screen-0.15.3.nf'

// DROP
include DROP from './NextflowModules/drop-1.2.4.nf'

def analysis_id = params.outdir.split('/')[-1]
def fastq_files = extractFastqPairFromDir(params.fastq_path)

workflow {
    // Use bam files if bam is set true
    if (!params.bam) {
        // FASTQC
        FastQC(fastq_files)

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
        Flagstat_raw( AlignReads.out.bam_file.join(Index.out).map {sample_id, rg_id, bam, bai ->
                                           [sample_id, bam, bai] } )
        bam_files = AlignReads.out.bam_file.map {sample_id, rg_id, bam ->
                                           [sample_id, bam] }
        bam_files = bam_files.join(Index.out)
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
    RSeQC_TIN(bam_files, genome_bed.collect())
    LCExtrap(bam_files)

    // Qualimap
    Qualimap(bam_files, genome_gtf)

    // Fastq screen - organism composition
    FastqScreen(fastq_files, params.fastq_screen_config)

    // MultiQC
    if (!params.bam) {
        qc_files = Channel.empty().mix(FastQC.out, TrimGalore.out.trimming_report, TrimGalore.out.fastqc_report, SortMeRNA.out.qc_report, star_logs, flagstat_logs, RSeQC.out, RSeQC_TIN.out, LCExtrap.out, Qualimap.out, PICARD_CollectMultipleMetrics.out, PICARD_EstimateLibraryComplexity.out, FastqScreen.out).collect()
    }
    else{
        qc_files = Channel.empty().mix(PICARD_CollectMultipleMetrics.out, PICARD_EstimateLibraryComplexity.out, RSeQC.out, RSeQC_TIN.out, LCExtrap.out, Qualimap.out, FastqScreen.out).collect()
    }
    MultiQC_post(analysis_id, qc_files)

    // WES data

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
            attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html"
        )

    } else {
        def subject = "RNAseq Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
