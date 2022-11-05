process CollectRnaSeqMetrics {
    label 'picard_2_27_5'
    container = 'broadinstitute/picard:2.27.5'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:

        path(bam_file)

    output:

      path("*.{RNA_Metrics,pdf}", emit: final_results_CollectRnaSeqMetrics)
      
    script:

	"""
	picard CollectRnaSeqMetrics \\
            --INPUT ${params.bam_infile}  \\
            --OUTPUT final_results_CollectRnaSeqMetrics \\
            --REF_FLAT ${params.REF_FLAT_infile} \\
            --STRAND NONE \\
            --CHART_OUTPUT chart.pdf
        """
}
