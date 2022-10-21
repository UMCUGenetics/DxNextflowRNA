process ASEReadCounter {
    tag {"ASEReadCounter"}
    label 'ASEReadCounter'
    container = 'broadinstitute/gatk:4.2.1.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(genome)
       	path(vcf_file)
       	path(bam_file)

    output:

        path("${vcf_file}")
        //path(vcf_file, emit: GATK_counts)
        //path("${vcf_infile}.txt")
        //path(".table", emit: final_results_ASEReadCounter_counts)

    script:

        """
        gatk ASEReadCounter \\
            -R ${params.genome}  \\
            -V ${params.vcf_infile}  \\
            -O final_results_ASEReadCounter_counts \\
            -I ${params.bam_infile}
        """
}


