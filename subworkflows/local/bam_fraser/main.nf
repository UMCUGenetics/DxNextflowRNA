include { FRASER            } from "../../../modules/local/fraser/main"
include { samplesheetToList } from 'plugin/nf-schema'

workflow BAM_FRASER {
    take:
    ch_bam_bai

    main:

    ch_refset = Channel.fromList(
        samplesheetToList(
            file(params.fraser_samplesheet),
            "${projectDir}/assets/refset_schema.json"))
        .map{ meta, basename ->
            junction_counts    = file("${params.fraser_reference_base}/${basename}_junction_counts.tsv.gz", checkIfExists: true)
            splice_site_counts = file("${params.fraser_reference_base}/${basename}_splice_counts.tsv.gz", checkIfExists: true)

            key = [ "${meta.treatment}_${meta.Set}"]

            [key, junction_counts, splice_site_counts]
        }
        .groupTuple() // group on treatment & set together
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
