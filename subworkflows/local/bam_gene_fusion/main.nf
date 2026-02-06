include { STARFUSION_DETECT } from '../../../modules/nf-core/starfusion/detect/main'


workflow BAM_GENE_FUSION {
    
    take:
    ch_star_junctions
    ch_starfusion_ref

    main:
    
    STARFUSION_DETECT(
        ch_star_junctions.map{ meta, junc -> [ meta, [], junc ]},
        ch_starfusion_ref
    )

    emit:
    starfusion_fusions = STARFUSION_DETECT.out.fusions
    starfusion_abridged = STARFUSION_DETECT.out.abridged
    starfusion_coding_effect = STARFUSION_DETECT.out.coding_effect ?: [:]
}
