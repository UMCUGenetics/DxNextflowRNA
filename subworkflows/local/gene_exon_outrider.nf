
include { OUTRIDER as OUTRIDER_EXON        } from '../../modules/local/outrider/main'
include { OUTRIDER as OUTRIDER_GENE        } from '../../modules/local/outrider/main'

workflow GENE_EXON_OUTRIDER {
    take:
    gene_counts // gene expression counts
    exon_counts // exon expression counts
    ch_gtf

    main:


    // ---------- helpers ----------
    // build channel with reference panels for each level ('gene'|'exon') and condition ('chx'|'cntrl')
    def get_refset = { String level, String cond ->
        parents = params."outrider_${level}_${cond}_parents"
        index   = params."outrider_${level}_${cond}_index"

        Channel.fromPath("${parents}/*")
            .collect()
            .map { ref -> [[id: "${cond}_parents"], ref] }
            .concat(
                Channel.fromPath("${index}/*")
                    .collect()
                    .map { ref -> [[id: "${cond}_index"], ref] }
            )
    }

    // Split CHX/control, group counts per treatment, combine them with ref panels
    def outriderInput = { countsCh, refsChx, refsCntrl ->
        def split = countsCh
            .map { meta, counts -> counts }
            .branch { v ->

                def s = v.toString().toLowerCase()
                def isChx   = s.contains('chx')
                def isCntrl = s.contains('control') || s.contains('cntrl') || s.contains('ctrl')

                chx  : isChx
                cntrl: isCntrl
            }

        def chx_in = split.chx
            .collect()
            .map { counts -> [counts] }
            .combine(refsChx)
            .map { counts, meta, ref -> [meta, counts, ref] }

        def cntrl_in = split.cntrl
            .collect()
            .map { counts -> [counts] }
            .combine(refsCntrl)
            .map { counts, meta, ref -> [meta, counts, ref] }

        chx_in.concat(cntrl_in)
    }


    //
    // OUTRIDER gene workflow
    //
    gene_refs_chx   = get_refset('gene','chx')
    gene_refs_cntrl = get_refset('gene','cntrl')
    gene_input      = outriderInput(gene_counts, gene_refs_chx, gene_refs_cntrl)

    OUTRIDER_GENE(gene_input, ch_gtf)

    //
    // EXON
    //
    exon_refs_chx   = get_refset('exon','chx')
    exon_refs_cntrl = get_refset('exon','cntrl')
    exon_input      = outriderInput(exon_counts, exon_refs_chx, exon_refs_cntrl)

    OUTRIDER_EXON(exon_input, ch_gtf)


    ch_versions = OUTRIDER_GENE.out.versions //OUTRIDER_EXON will be the same version

    emit:
    gene_tsv_full   = OUTRIDER_GENE.out.tsv_full // channel: [ val(meta), path(*_outrider_results.tsv) ]
    gene_tsv_signif = OUTRIDER_GENE.out.tsv_signif // channel: [ val(meta), path(*_outrider_results.tsv) ]
    exon_tsv_full   = OUTRIDER_EXON.out.tsv_full // channel: [ val(meta), path(*_outrider_results.tsv) ]
    exon_tsv_signif = OUTRIDER_EXON.out.tsv_signif // channel: [ val(meta), path(*_outrider_results.tsv) ]
    gene_multiqc    = OUTRIDER_GENE.out.multiqc_tsv
    exon_multiqc    = OUTRIDER_EXON.out.multiqc_tsv
    versions        = ch_versions

}
