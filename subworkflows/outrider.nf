#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { OUTRIDER } from '../NextflowModules/Outrider/1.20.0/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outrider (sub)workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow outrider {
    take:
	featurecounts

    main:
        ch_outrider_ref_gene = params.refgene.contains(",") ? Channel.fromPath(params.refgene?.split(',') as List).view() : Channel.fromPath("$params.refgene").view()
        ch_outrider_ref_exon = params.refexon.contains(",") ? Channel.fromPath(params.refexon?.split(',') as List).view() : Channel.fromPath("$params.refexon").view()

        OUTRIDER(
	    featurecounts,
            ch_outrider_ref_gene.collect(),
            ch_outrider_ref_exon.collect()
        )
}


workflow outrider_entry {
    featurecounts_gene = Channel.fromPath("$params.outdir/feature_counts/*gene*.txt").map { counts -> [[id:counts.getSimpleName()], counts, "gene"] }
    featurecounts_exon = Channel.fromPath("$params.outdir/feature_counts/*exon*.txt").map { counts -> [[id:counts.getSimpleName()], counts, "exon"] }

    outrider(featurecounts_gene.concat(featurecounts_exon))
}