#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// fastq modules
include extractFastqPairFromDir from './NextflowModules/Utils/fastq.nf'

//QC module
include FastQC from './NextflowModules/FastQC/0.11.8/FastQC.nf' params(optional:'')
include MultiQC from './NextflowModules/MultiQC/1.10/MultiQC.nf' params(optional:'')

def fastq_files = extractFastqPairFromDir(params.fastq_path)
def analysis_id = params.outdir.split('/')[-1]

workflow {

    // FASTQC
    FastQC(fastq_files)

    // MultiQC
    MultiQC(analysis_id, FastQC.out)
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
        def subject = "WES Workflow Successful: ${analysis_id}"
        sendMail(
            to: params.email.trim(),
            subject: subject,
            body: email_html,
            attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html"
        )

    } else {
        def subject = "WES Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
