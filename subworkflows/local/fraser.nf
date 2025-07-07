include { FRASER } from "../../modules/local/fraser/main"

workflow RARE_SPLICING {
    take:
    ch_bam_bai
    refset_path


    main:
    ch_refset = Channel.fromFilePairs(
            refset_path,
            checkIfExists: true)
            { file -> file.name.replaceAll(/.bam|.bai$/,'') }
                .map{ meta, bam_index -> [bam_index[0], bam_index[1]] }
                .reduce( [[], []] ) { acc, pair ->
                    acc[0] << pair[0]      // BAM files
                    acc[1] << pair[1]      // BAI files
                    acc                    // final output
            }

    FRASER(
        ch_bam_bai,
        ch_refset
    )

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(FRASER.out.versions) //OUTRIDER_EXON will be the same version

    emit:
    tsv      = FRASER.out.tsv
    heatmap  = FRASER.out.heatmap
    versions = ch_versions

}
