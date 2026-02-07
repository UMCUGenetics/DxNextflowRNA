include { STARFUSION_DETECT } from '../../../modules/nf-core/starfusion/detect/main'
include { ARRIBA_ARRIBA } from '../../../modules/nf-core/arriba/arriba/main' 


workflow BAM_GENE_FUSION {
    
    take:
    ch_star_junctions
    ch_starfusion_ref
    ch_bam
    ch_genome_fasta
    ch_genome_gtf
        

    main:
    
    STARFUSION_DETECT(
        ch_star_junctions.map{ meta, junc -> [ meta, [], junc ]},
        ch_starfusion_ref
    )

    ARRIBA_ARRIBA(
        ch_bam,
        ch_genome_fasta,
        ch_genome_gtf,
        [], //blacklist
        [], //known fusions
        [], //cytobands
        [] //protein domains
    )

    emit:
    starfusion_fusions       = STARFUSION_DETECT.out.fusions
    starfusion_abridged      = STARFUSION_DETECT.out.abridged
    starfusion_coding_effect = STARFUSION_DETECT.out.coding_effect ?: [:]
    ariba_fusions            = ARRIBA_ARRIBA.out.fusions
}
