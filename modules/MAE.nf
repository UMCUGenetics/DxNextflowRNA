process MAE {

    tag {"MonoAllelicExpression"}
    label 'MAE'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(table_file)

    output:
        path(MAE_output)
        path("*.{csv}", emit: MAE_output)

    script:
        """
        Rscript ${baseDir}/MAE.r ${table_file} ${MAE_output}
        """
}
