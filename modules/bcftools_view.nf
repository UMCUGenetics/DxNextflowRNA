process view_bcf_norm {

    tag {"BCFtools View_BCF_VCF ${vcf_file}"}
    label 'BCFtools_1_10_2'
    label 'BCFtools_1_10_2_View_BCF_VCF'
    container = 'quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    file(vcf_file)

    output:
    file("${vcf_file}.vcf")

    script:
    """
    bcftools view ${vcf_file} ${params.optional} > ${vcf_file}.vcf
    """
}
