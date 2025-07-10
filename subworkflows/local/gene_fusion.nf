include { STARFUSION                    } from '../../modules/local/star/fusion/main'

workflow GENE_FUSION {
    take:
    ch_junction
    ch_star_fusion_index

    main:
    STARFUSION(
        ch_junction,
        ch_star_fusion_index
    )

    emit:
    versions = STARFUSION.out.versions
}