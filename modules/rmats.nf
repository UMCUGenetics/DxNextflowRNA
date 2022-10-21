process RMATS {
    tag {"RMATS"}
    label 'RMATS_412'
    container = 'quay.io/biocontainers/rmats:4.1.2--py37haf75f70_1'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(s1)
        path(genome_gtf)
        path(STAR_genome_index)

    output:
        path("*", emit: final_results_rmats)


    script:
    """
    rmats.py \
      --gtf gencode.v39.annotation.gtf \\
      --bi hg38_genome \\
      --readLength 150 \\
      --nthread 10 \\
      --novelSS \\
      --statoff \\
      --mil 50 \\
      --mel 500 \\
      --s1 ${s1} \\
      --od final_results_rmats \\
      --tmp final_results_rmats

    """
 }
