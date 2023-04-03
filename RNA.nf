#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq module
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

// QC modules
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MultiQC as MultiQC_pre from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'')
include MultiQC as MultiQC_post from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'')

// Trimming module
include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params(optional:'--fastqc', singleEnd: params.singleEnd)

// Mapping modules
include GenomeGenerate from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf'
include AlignReads from './NextflowModules/STAR/2.7.3a/AlignReads.nf'
include Index from './NextflowModules/Sambamba/0.7.0/Index.nf'
include Flagstat as Flagstat_raw from './NextflowModules/Sambamba/0.7.0/Flagstat.nf'

// DROP
include DROP from './Drop-1.2.4.nf'

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

workflow {

    // FASTQC
    FastQC(fastq_files)

    // MultiQC
    MultiQC_pre(analysis_id, FastQC.out)

    // TrimGalore
    TrimGalore(fastq_files)
    trimmed_fastqs = TrimGalore.out.fastqs_trimmed

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
    AlignReads( fastq_files, star_index.collect(), genome_gtf.collect() )
    Index( AlignReads.out.bam_file.map {sample_id, rg_id, bam ->
                                       [sample_id, bam] })
    Flagstat_raw( AlignReads.out.bam_file.join(Index.out).map {sample_id, rg_id, bam, bai ->
                                       [sample_id, bam, bai] } )
    bam_sorted = AlignReads.out.bam_file.join(Index.out)
    star_logs = AlignReads.out.log.mix(AlignReads.out.final_log)
    flagstat_logs = Flagstat_raw.out

    // MultiQC
    qc_files = Channel.empty().mix(FastQC.out, TrimGalore.out.trimming_report, TrimGalore.out.fastqc_report, star_logs, flagstat_logs).collect()
    MultiQC_post(analysis_id, qc_files)

    // DROP
    DROP()
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
