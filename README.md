# UMCUGenetics-DxNextflowRNA
To run the workflow please use the following command:

`<path to nextflow bash script> <path to fastq directory> <path to WES directory> <path to output directory> <email adress>`

This project contains a Nextflow workflow for analyzing RNA-seq data to detect deviate RNA.

When the workflow is run with `bam=true` in `nextflow.config` then there should be bam files in the fastq directory as well as fastq files. 
 
## How to install DROP

```
wget https://github.com/nickhsmith/drop_demo_data/archive/refs/heads/main.zip
unzip main.zip
mv drop_demo_data-main/Data Data
rm -rf drop_demo_data-main
mkdir R_lib
singularity shell -B $TMPDIR:$TMPDIR -B /hpc:/hpc docker://quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0
export R_LIBS=R_lib:/usr/local/lib/R/library
drop demo
snakemake -n # dryrun
snakemake –cores 2 #run on test dataset
```

## How to run DROP on other data

Change the config.yaml inside the drop directory, focus on the 'groups' and make sure that these are identical to the 'groups' of the sample annotation file. 
Make sure that there are enough samples in your sample annotation, otherwise it will crash with statistics (splicing module arount 92%). An example of this is shown below:

```
Quitting from lines 73-77 (/scratch/212850/RtmpYjGn8c/file2375f71d0a0767/Summary.spin.Rmd) 
Error: BiocParallel errors
  1 remote errors, element index: 1
  0 unevaluated and other errors
  first remote error:
Error in .local(object, ...): Missing rawOtherCounts for type 'jaccard'. Please calculate PSIValues first. And then try again.
In addition: Warning message:
In .local(object, ...) : no hyperparameters were estimated for jaccard 
Please use `optimHyperParams` to compute them.

Execution halted
[Tue Apr 18 17:41:36 2023]
Error in rule AberrantSplicing_pipeline_FRASER_Summary_R:
    jobid: 32
    input: /hpc/diaggen/projects/RNAseq_Jade/dropOutput/processed_results/aberrant_splicing/datasets/savedObjects/fraser--v43/fds-object.RDS, /hpc/diaggen/projects/RNAseq_Jade/dropOutput/processed_results/aberrant_splicing/results/v43/fraser/fraser/results.tsv, Scripts/AberrantSplicing/pipeline/FRASER/Summary.R
    output: /hpc/diaggen/projects/RNAseq_Jade/dropOutput/html/AberrantSplicing/fraser--v43_summary.html, /hpc/diaggen/projects/RNAseq_Jade/dropOutput/html/AberrantSplicing/FRASER_results_fraser--v43.tsv
    log: /hpc/diaggen/projects/RNAseq_Jade/drop/.drop/tmp/AS/fraser--v43/FRASER_summary.Rds (check log file(s) for error message)

RuleException:
CalledProcessError in line 248 of /scratch/212850/tmpjor6rtv4:
Command 'set -euo pipefail;  Rscript --vanilla /hpc/diaggen/projects/RNAseq_Jade/drop/.snakemake/scripts/tmpxfu889z3.wBRender.R' returned non-zero exit status 1.
  File "/scratch/212850/tmpjor6rtv4", line 248, in __rule_AberrantSplicing_pipeline_FRASER_Summary_R
  File "/usr/local/lib/python3.11/concurrent/futures/thread.py", line 58, in run
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2023-04-18T104052.676729.snakemake.log
``` 

When you have adapted the config.yaml, create a file sample_annotation.tsv, make sure that the path is stated in the config.yaml. 
The sample annotation file should look like this:

```
RNA_ID	RNA_BAM_FILE	DROP_GROUP	PAIRED_END	COUNT_MODE	COUNT_OVERLAPS	DNA_VCF_FILE	DNA_ID	STRAND	TISSUE	SEX
ERR188023	/hpc/diaggen/projects/RNAseq_Jade/data/geuvadis/RNAseq_bams/ERR188023_null_0_Aligned.sortedByCoord.out.bam	blood	TRUE	IntersectionStrict	TRUE			no	Blood	female
ERR188024	/hpc/diaggen/projects/RNAseq_Jade/data/geuvadis/RNAseq_bams/ERR188024_null_0_Aligned.sortedByCoord.out.bam	blood	TRUE	IntersectionStrict	TRUE			no	Blood	female
2023-803CHX	/hpc/diaggen/projects/RNAseq_Jade/patienten_data/RNAseq_bams/2023-803CHX.bam	blood,mae	TRUE	IntersectionStrict	TRUE	/hpc/diaggen/projects/RNAseq_Jade/patienten_data/vcf/U189155CF2021D20023_patienten_data_RNA_WES.vcf.gz	U189155CF2021D20023	no	Fibroblasts	female
2023-803Cntrl	/hpc/diaggen/projects/RNAseq_Jade/patienten_data/RNAseq_bams/2023-803Cntrl.bam	blood,mae	TRUE	IntersectionStrict	TRUE	/hpc/diaggen/projects/RNAseq_Jade/patienten_data/vcf/U189155CF2021D20023_patienten_data_RNA_WES.vcf.gz	U189155CF2021D20023	no	Fibroblasts	female
...
```

This file has to contain at least 10 samples before DROP will run, but will then also probably crash at the statistics. 

## To do:

- Automatically creating a sample annotation so that DROP can be run inside the workflow

- Adding STARfusion

- Making every tool flexible on time and memory
