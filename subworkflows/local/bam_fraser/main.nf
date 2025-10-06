include { FRASER             } from "../../../modules/local/fraser/run/main"
include { FRASER_BUILDCOUNTS } from "../../../modules/local/fraser/counts/main"
include { FRASER_CALC_COUNTS } from "../../../modules/local/fraser/calc_counts/main"
include { samplesheetToList  } from 'plugin/nf-schema'

workflow BAM_FRASER {

    take:
    ch_bam_bai

    main:

    ch_refset = Channel.fromList(
        samplesheetToList(
            file(params.fraser_samplesheet),
            "${projectDir}/assets/refset_schema.json"))
        .map{ meta, basename ->
            def junction_counts = file("${params.fraser_reference_base}/${basename}_junction_counts.tsv.gz", checkIfExists: true)
            def key =  "${meta.treatment}_${meta.Set}"

            [key, junction_counts]
        }
        .groupTuple() // group on treatment & set together
        .map { group_id, junction_counts ->
            def meta = [id: group_id]
            [meta, junction_counts]
        }

    FRASER_CALC_COUNTS(
        ch_bam_bai
    )


    FRASER_BUILDCOUNTS(
        ch_refset,
        FRASER_CALC_COUNTS.out.junctionCounts.collect()
    )

    FRASER(
        FRASER_BUILDCOUNTS.out.countTable
            .join(FRASER_BUILDCOUNTS.out.spliceCounts)
            .join(FRASER_BUILDCOUNTS.out.junctionCounts)
    )

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(FRASER.out.versions)

    emit:
    tsv      = FRASER.out.tsv
}
