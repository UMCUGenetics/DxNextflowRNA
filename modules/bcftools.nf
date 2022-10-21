process norm_bcf_norm {

    tag {"BCFtools Norm_BCF_VCF ${vcf_file}"}
    label 'BCFtools_1_10_2'
    label 'BCFtools_1_10_2_Norm_BCF_VCF'
    container = 'quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(vcf_file)

    output:
        path("${vcf_file}.vcf")

    script:
        """
        bcftools norm ${params.optional} ${vcf_file} > ${vcf_file}.vcf
        """
}
