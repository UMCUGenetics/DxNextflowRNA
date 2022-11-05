process RMATS {

    label 'rmats_4.1.2'
    label 'rmats_4.1.2_RMATS_AS'
    container = 'quay.io/biocontainers/rmats:4.1.2--py37haf75f70_1'
    shell = ['/bin/bash', '-euo', 'pipefail']


    input:
    path(STAR_genome_index)
    path(genome_gtf)
    path(s1)
    path(s2)

    output:
    path("*.txt", emit: final_results_rmats)


    script:
    """
  rmats.py \
      --gtf gencode.v39.annotation.gtf \\
      --bi hg38_genome \\
      --s2 s2.txt \\
      --readLength 150 \\
      --nthread 10 \\
      --novelSS \\
      --statoff \\
      --mil 50 \\
      --mel 500 \\
      --s1 s1.txt \\
      --od final_results_rmats \\
      --tmp final_results_rmats

    """
 }
