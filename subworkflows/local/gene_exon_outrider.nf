
include { OUTRIDER as OUTRIDER_EXON        } from '../../modules/local/outrider/main'
include { OUTRIDER as OUTRIDER_GENE        } from '../../modules/local/outrider/main'

workflow GENE_EXON_OUTRIDER {
    take:
    gene_counts // gene expression counts
    exon_counts // exon expression counts
    ch_gtf

    main:
    //
    // OUTRIDER gene workflow
    //
    ch_outrider_gene_chx_set1 = Channel.fromPath("${params.outrider_gene_chx_parents}/*")
        .collect()
        .map{ ref -> [[id:"chx_parents"], ref] }
    ch_outrider_gene_chx_set2 = Channel.fromPath("${params.outrider_gene_chx_index}/*")
        .collect()
        .map{ ref -> [[id:"chx_index"], ref] }
    ch_outrider_gene_cntrl_set1 = Channel.fromPath("${params.outrider_gene_cntrl_parents}/*")
        .collect()
        .map{ ref -> [[id:"cntrl_parents"], ref] }
    ch_outrider_gene_cntrl_set2 = Channel.fromPath("${params.outrider_gene_cntrl_index}/*")
        .collect()
        .map{ ref -> [[id:"cntrl_index"], ref] }

    ch_outrider_gene_chx = ch_outrider_gene_chx_set1
        .concat(ch_outrider_gene_chx_set2)

    ch_outrider_gene_cntrl = ch_outrider_gene_cntrl_set1
        .concat(ch_outrider_gene_cntrl_set2)

    ch_feature_counts_genes = gene_counts
            .map { meta, counts -> counts }
            .branch { v ->
                chx: v.toString().contains('CHX')
                cntrl: v.toString().contains("cntrl")
            }

    ch_outrider_gene_chx_input = ch_feature_counts_genes.chx
        .collect()
        .map{ counts -> [counts] }
        .combine(ch_outrider_gene_chx)
        .map{ counts, meta, ref -> [meta, counts, ref] }

    ch_outrider_gene_cntrl_input = ch_feature_counts_genes.cntrl
        .collect()
        .map{ counts -> [counts] }
        .combine(ch_outrider_gene_cntrl)
        .map{ counts, meta, ref -> [meta, counts, ref] }

    ch_outrider_gene_input = ch_outrider_gene_chx_input
        .concat(ch_outrider_gene_cntrl_input)

    OUTRIDER_GENE(
        ch_outrider_gene_input,
        ch_gtf
    )

    //
    // OUTRIDER exon workflow
    //
    ch_feature_counts_exons = exon_counts
            .map { meta, counts -> counts }
            .branch { v ->
                chx: v.toString().contains('CHX')
                cntrl: v.toString().contains("cntrl")
            }

    ch_outrider_exon_chx_set1 = Channel.fromPath("${params.outrider_exon_chx_parents}/*")
        .collect()
        .map{ ref -> [[id:"chx_parents"], ref] }
    ch_outrider_exon_chx_set2 = Channel.fromPath("${params.outrider_exon_chx_index}/*")
        .collect()
        .map{ ref -> [[id:"chx_index"], ref] }
    ch_outrider_exon_cntrl_set1 = Channel.fromPath("${params.outrider_exon_cntrl_parents}/*")
        .collect()
        .map{ ref -> [[id:"cntrl_parents"], ref] }
    ch_outrider_exon_cntrl_set2 = Channel.fromPath("${params.outrider_exon_cntrl_index}/*")
        .collect()
        .map{ ref -> [[id:"cntrl_index"], ref] }

    ch_outrider_exon_chx = ch_outrider_exon_chx_set1
        .concat(ch_outrider_exon_chx_set2)

    ch_outrider_exon_cntrl = ch_outrider_exon_cntrl_set1
        .concat(ch_outrider_exon_cntrl_set2)

    ch_outrider_exon_chx_input = ch_feature_counts_exons.chx
        .collect()
        .map{ counts -> [counts] }
        .combine(ch_outrider_exon_chx)
        .map{ counts, meta, ref -> [meta, counts, ref] }

    ch_outrider_exon_cntrl_input = ch_feature_counts_exons.cntrl
        .collect()
        .map{ counts -> [counts] }
        .combine(ch_outrider_exon_cntrl)
        .map{ counts, meta, ref -> [meta, counts, ref] }

    ch_outrider_exon_input = ch_outrider_exon_chx_input
        .concat(ch_outrider_exon_cntrl_input)

    OUTRIDER_EXON(
        ch_outrider_exon_input,
        ch_gtf
    )

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(OUTRIDER_GENE.out.versions) //OUTRIDER_EXON will be the same version


    emit:
    gene_tsv_full   = OUTRIDER_GENE.out.tsv_full // channel: [ val(meta), path(*_outrider_results.tsv) ]
    gene_tsv_signif = OUTRIDER_GENE.out.tsv_signif // channel: [ val(meta), path(*_outrider_results.tsv) ]
    gene_volcano    = OUTRIDER_GENE.out.volcano_plots // channel: [ val(meta), path(*_volcano_plots/) ]
    exon_tsv_full   = OUTRIDER_EXON.out.tsv_full // channel: [ val(meta), path(*_outrider_results.tsv) ]
    exon_tsv_signif = OUTRIDER_EXON.out.tsv_signif // channel: [ val(meta), path(*_outrider_results.tsv) ]
    exon_volcano    = OUTRIDER_EXON.out.volcano_plots // channel: [ val(meta), path(*_volcano_plots/) ]
    gene_multiqc    = OUTRIDER_GENE.out.multiqc_tsv
    exon_multiqc    = OUTRIDER_EXON.out.multiqc_tsv
    versions        = ch_versions

}
