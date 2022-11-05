process bgzip {

    tag {"Samtools BGZIP ${vcf_file}"}
    label 'Samtools_1_10'
    label 'Samtools_1_10_BGZIP'
    container = 'quay.io/biocontainers/samtools:1.10--h9402c20_2'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    file(vcf_file)

    output:
    file("${vcf_file}.vcf.gz")

    script:
    """
    bgzip ${params.optional} ${vcf_file} > ${vcf_file}.vcf.gz
    """
}

