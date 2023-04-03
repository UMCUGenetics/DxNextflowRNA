process DROP {
    tag {"DROP ${sample_id} - ${rg_id}"}
    label 'DROP_1_2_4'
    container = 'docker://quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple (sample_id, rg_id )

    script:
        """
        export R_LIBS=${params.drop_path}/R_lib:/usr/local/lib/R/library
        drop demo
        snakemake -n > test_drop.txt
        """
}
