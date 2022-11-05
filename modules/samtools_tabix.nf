process tabix {

    tag {"Samtools ${gz_vcf_file}"}
    label 'Samtools_1_10'
    container = 'quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    file(gz_vcf_file)


    script:
    """
    tabix ${gz_vcf_file}
    """
}
