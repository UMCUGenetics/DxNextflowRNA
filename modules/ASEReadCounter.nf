process ASEReadCounter {
//    tag {"GATK GenotypeGVCF ${sample_id} - ${interval_file.simpleName}"}
//    label 'GATK_4_2_1_0'
//    label 'GATK_4_2_1_0_GenotypeGVCF'
//    container = 'broadinstitute/gatk:4.2.1.0'
//    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
//        tuple(sample_id, path(gvcf_files), path(gvcf_idx_files), path(interval_file))
       	path(vcf_file)
       	path(bam_file)
        path(genome)
        

    output:

        path("${vcf_file}")
//        path(vcf_file, emit: GATK_counts)
//        path("${vcf_infile}.txt")
//        path(".table", emit: final_results_ASEReadCounter_counts)

    script:

        """
        gatk ASEReadCounter \\
            -R ${params.genome}  \\
            -V ${params.vcf_infile}  \\
            -O final_results_ASEReadCounter_counts \\
            -I ${params.bam_infile}
        """
}


// gatk ASEReadCounter  -R /hpc/diaggen/users/Behzad/bam/STAR/hg38_genome/GRCh38.p13.genome.fa  -I PMABM000INO_PMCRZ304SHHAligned.out.sorted_2.bam  -V PMABM000INO_PMCRZ629LVV_RNA-Seq.vcf.gz  -O output3.table
