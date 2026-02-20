include { STARFUSION_DETECT     } from '../../../modules/nf-core/starfusion/detect/main'
include { AGAT_CONVERTSPGFF2TSV } from '../../../modules/nf-core/agat/convertspgff2tsv/main'
include { ARRIBA_ARRIBA         } from '../../../modules/nf-core/arriba/arriba/main'
include { FUSIONINSPECTOR       } from '../../../modules/nf-core/fusioninspector/main'
include { FUSIONREPORT_DETECT   } from '../../../modules/nf-core/fusionreport/detect/main'
include { FUSIONREPORT_DOWNLOAD } from '../../../modules/nf-core/fusionreport/download/main'
include { VCF_COLLECT           } from '../../../modules/local/vcf_collect/main'
include { HGNC_DOWNLOAD         } from '../../../modules/local/hgnc_download/main'

workflow BAM_GENE_FUSION {
    take:
    ch_star_junctions //: Channel<Tuple<Map, Path>>
    ch_starfusion_ref
    ch_bam
    ch_fastq
    ch_genome_fasta
    ch_genome_gtf
    arriba_blacklist
    arriba_known_fusions
    arriba_cytobands
    arriba_protein
        

    main:
    
    STARFUSION_DETECT(
        ch_star_junctions.map{ meta, junc -> [ meta, [], junc ]},
        ch_starfusion_ref
    )

    ARRIBA_ARRIBA(
        ch_bam.map{ meta, bam, _bai -> [meta, bam] },
        ch_genome_fasta.map{ meta, fasta, _fai -> [meta, fasta] },
        ch_genome_gtf,
        arriba_blacklist, 
        arriba_known_fusions, 
        arriba_cytobands, 
        arriba_protein 
    )

    def ch_fusions = ARRIBA_ARRIBA.out.fusions
        .join(STARFUSION_DETECT.out.fusions, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, arriba, starfusion -> [meta, arriba, starfusion, []] }

    FUSIONREPORT_DOWNLOAD()

    ch_fusionreport_db = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref.collect()

    FUSIONREPORT_DETECT(
        ch_fusions,
        ch_fusionreport_db,
        params.fusion_tools_cutoff
    )

    ch_reads_fusions = ch_fastq
        .map{ meta, reads -> [meta, reads] }
        .join(FUSIONREPORT_DETECT.out.fusion_list
              .map{meta, fusion_list -> [meta,[fusion_list]]})
    
    FUSIONINSPECTOR(
        ch_reads_fusions,
        ch_starfusion_ref.map{fusion_ref ->
            [[id: "starfusion_index"], fusion_ref]}
    )

    HGNC_DOWNLOAD()

    AGAT_CONVERTSPGFF2TSV(
        FUSIONINSPECTOR.out.out_gtf
            .filter { _meta, file -> file.exists() && file.size() > 0 }
    )
    
    VCF_COLLECT(
        FUSIONINSPECTOR.out.abridged_tsv
            .filter{ _meta, file -> file.exists() && file.size() > 0 }
            .join(AGAT_CONVERTSPGFF2TSV.out.tsv)
            .join(FUSIONREPORT_DETECT.out.report)
            .join(FUSIONREPORT_DETECT.out.csv),
        HGNC_DOWNLOAD.out.hgnc_ref
            .map{ hgnc_ref -> [[id: hgnc_ref.getSimpleName()], hgnc_ref]},
        HGNC_DOWNLOAD.out.hgnc_date
            .map{ hgnc_date -> [[id: hgnc_date.getSimpleName()], hgnc_date]}  
    )

    
    emit:
    starfusion_fusions       = STARFUSION_DETECT.out.fusions
    starfusion_abridged      = STARFUSION_DETECT.out.abridged
    starfusion_coding_effect = STARFUSION_DETECT.out.coding_effect ?: [:]
    ariba_fusions            = ARRIBA_ARRIBA.out.fusions
}
