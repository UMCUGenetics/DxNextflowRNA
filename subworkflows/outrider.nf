#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { OUTRIDER as OUTRIDER_GENE } from '../NextflowModules/Outrider/1.20.0/main'
include { OUTRIDER as OUTRIDER_EXON } from '../NextflowModules/Outrider/1.20.0/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outrider subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow outrider {
    take:
	    featurecounts_gene
	    featurecounts_exon

    main:
    ch_gtf = Channel.fromPath(params.gtf).map { gtf -> [gtf.getSimpleName(), gtf] }.first()
	if (params.refgene != null && params.inputgene != null){
	    ch_outrider_ref_gene = params.refgene.contains(",") ? Channel.fromPath(params.refgene?.split(',') as List) : Channel.fromPath("$params.refgene")

            OUTRIDER_GENE(
	            featurecounts_gene,
                ch_outrider_ref_gene.collect(),
                ch_gtf
            )
	}

	if (params.refexon != null && params.inputexon != null){
            ch_outrider_ref_exon = params.refexon.contains(",") ? Channel.fromPath(params.refexon?.split(',') as List) : Channel.fromPath("$params.refexon")
            OUTRIDER_EXON(
                featurecounts_exon,
                ch_outrider_ref_exon.collect(),
                ch_gtf
            )
	}
}


workflow outrider_entry {
    featurecounts_gene = Channel.fromPath("$params.inputgene").map { counts -> [[id:counts.getSimpleName()], counts] }
    featurecounts_exon = Channel.fromPath("$params.inputexon").map { counts -> [[id:counts.getSimpleName()], counts] }
    outrider(featurecounts_gene, featurecounts_exon)
}
