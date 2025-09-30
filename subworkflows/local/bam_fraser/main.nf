include { FRASER            } from "../../../modules/local/fraser/main"
include { samplesheetToList } from 'plugin/nf-schema'

workflow BAM_FRASER {
    take:
    ch_bam_bai


    main:

    ch_refset = Channel.fromList(
        samplesheetToList(file(params.fraser_samplesheet), "${projectDir}/assets/refset_schema.json"))
        .map{ meta, bam ->
            bai = file("${bam}.bai", checkIfExists: true)
            def key = [ "${meta.treatment}_${meta.Set}"]

            [key, file(bam, checkIfExists: true), bai]
        }
        .groupTuple()
        .map { group_id, bams, bais ->
            def meta = [id: group_id]

            [meta, bams, bais]
        }

    FRASER(
        ch_bam_bai,
        ch_refset
    )

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(FRASER.out.versions)

    emit:
    tsv      = FRASER.out.tsv
}
