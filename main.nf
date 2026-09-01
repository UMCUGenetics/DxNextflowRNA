#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    umcugenetics/dxnextflowrna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/umcugenetics/dxnextflowrna
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_VARIANT_CALLING     } from './subworkflows/local/bam_variant_calling/main'
include { DXNEXTFLOWRNA           } from './workflows/dxnextflowrna'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_umcugenetics_dxnextflowrna_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_umcugenetics_dxnextflowrna_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
    )

    //
    // WORKFLOW: Run main workflow
    //
    UMCUGENETICS_DXNEXTFLOWRNA()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        UMCUGENETICS_DXNEXTFLOWRNA.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow UMCUGENETICS_DXNEXTFLOWRNA {
    main:

    //
    // WORKFLOW: Run pipeline
    //
    DXNEXTFLOWRNA()

    emit:
    multiqc_report = DXNEXTFLOWRNA.out.multiqc_report // channel: /path/to/multiqc_report.html
}

workflow DXNEXTFLOWRNA_VARIANT_CALLING {
    main:

    PREPARE_REFERENCES()

    ch_bam_bai = Channel.fromFilePairs(
        "${params.input}/*.{bam,bai}",
        checkIfExists: true
    ) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
        .map{ meta, bam_index -> [['id': meta], bam_index[0], bam_index[1]] }

    BAM_VARIANT_CALLING(
        ch_bam_bai,
        PREPARE_REFERENCES.out.ch_fasta_fai,
        PREPARE_REFERENCES.out.dict,
        PREPARE_REFERENCES.out.interval_list_split,
        PREPARE_REFERENCES.out.ch_dbsnp,
        PREPARE_REFERENCES.out.ch_dbsnp_tbi
    )

}
