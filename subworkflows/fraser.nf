#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FRASER } from '../NextflowModules/Fraser/1.99.3/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Fraser subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow fraser {
    take:
            input_bam

    main:
	if (params.refbam != null && params.inputbam != null){
            ch_ref_bam = params.refbam.contains(",") ? Channel.fromPath(params.refbam?.split(',') as List) : Channel.fromPath("$params.refbam")

            FRASER(
                input_bam,
                ch_ref_bam.collect(),
            )
	}

}


workflow fraser_entry {
    ch_input_bam = Channel.fromPath("$params.inputbam").map { inputbam -> [[id:inputbam.getSimpleName()], inputbam] }
    fraser(ch_input_bam)
}
