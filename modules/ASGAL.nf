process ASGAL {

    label 'ASGAL_v1.1.7'
    label 'ASGAL_v1.1.7_Novel_AS'
    container = 'algolab/asgal:v1.1.7'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    path genome
    path transcriptome
    path gtf
    tuple val(pair_id),rg_id, path(reads)


    output:

    path("*.{mem,sam,csv}", emit: final_ASGAL)


    script:
    def barcode = rg_id.split('_')[1]

    """
    python3 /galig/asgal \\
        ${params.optional} \\
        -g ${genome}\\
        -t ${transcriptome}  \\
        -a ${gtf}  \\
        -s ${reads[0]} \\
        -s2 ${reads[1]} \\
        -@ 1 \\
        -o final_ASGAL
    """
}
