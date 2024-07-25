[![GitHub Actions CI Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20CI/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20linting/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# DxNextflowRNA
UMCU Genetics Nextflow RNA workflow

## Usage

```bash
export NXF_JAVA_HOME='/hpc/diaggen/projects/woc/rna/tools/jdk-18.0.2.1'
/hpc/diaggen/projects/woc/rna/tools/nextflow run DxNextflowRNA --input /hpc/diaggen/projects/woc/rna/testdata/fastq --outdir ...
```

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  DxNextflowRNA for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


# Instructions

## Workflows
Instruction to use the pipeline
HPC locations for reference sets are provided at the SOP

1.	To run the *full workflow*, reference sets for diagnostic pipeline use, gene and exon should be used

Start workflow with:

```bash
bash run_nextflow_RNAseq.sh {input fastq directory} {output directory} {email} --refgene { dir_refset_genelevel} --refexon {dir_refset_genelevel}
```


Fraser can also be added to the full workflow:

DxNextflowRNA, main.nf, comment out: 
```bash
//fraser(fastq_to_bam.out)
```
use:
```bash
bash run_nextflow_RNAseq.sh {input fastq directory} {output directory} {email} --refgene { dir_refset_genelevel} --refexon {dir_refset_genelevel} --refbam{ dir_bam_bambai_referenceset}
```

2.	**subworkflow fastq to bam**. Bam files are collected at {output directory}/bam_files/

```bash
bash run_nextflow_RNAseq.sh {input fastq directory} {output directory} {email} -entry fastq_to_bam
```

3.	**subworkflow featurecounts**: {input bam_files directory} is used for bam files input directory. Bam files (can also be) copied to {output directory}/bam_files with symlink (ln -s). 

Featurecounts files are collected at {output directory}/feature_counts/

Start workflow: 

```bash
bash run_nextflow_RNAseq.sh {bamfiles directory} {output directory} {email} -entry featurecounts
```

4.	**subworkflow outrider**: {input fastq directory} is not used (. can be provided), featurecounts files for sample of interest and reference set need to be provided. 

Outrider output files collected at {output directory}/outrider/

Start workflow:

```bash
bash run_nextflow_RNAseq.sh {input fastq directory} {output directory} {email} -entry outrider_entry --inputgene {featurecounts_file_gene_sampleofinterest} --refgene {dir_refset_genelevel} --inputexon {featurecounts_file_exon_sampleofinterest} --refexon {outrider/refset_exonlevel}
```


The outrider subworkflow generates appropriate files for Jupyter webapp automatically

5.	**subworkflow fraser**: {input fastq directory} is not used (. Can be provided), bam/bam.bai directories for sample of interest and reference set need to be provided. 

```bash
bash run_nextflow_RNAseq.sh {input fastq directory} {output directory} {email} -entry fraser_entry --inputbam {dir_bam_bambai_sampleofinterest} --refbam { dir_bam_bambai_referenceset}
```


## Convert data for EMC APP

1.	Countdata, metadata and outrider result files are collected in {output_directory}/outrider. 


umcu_rnaseq_metadata.csv is identical for gene and exon level. Chx metadata file should be added to cntrl(untreated) metadata file: 

```bash
cat {metadata cntrl file} < (tail -n+2 {metadata chx file} )  > umcu_rnaseq_metadata.csv
```

Add this file to metadata directory of EMC App: RNASEQ-VOILA/metadata


cntrl and chx countdata files pasted to 1 file for gene and one for exon. 

Use script/merge_countdata.sh, options -c {cntrl file in} -h {chx file in} -o {gene or exon}

A file umcu_rnaseq_{gene/exon}_counts.tsv is created. 

Copy this file to RNASEQ-VOILA/countdata/

2.	Add gene associated phenotypes to outrider result files: umcu_rnaseq_fib_{untreated/chx}_res_outrider_{gene/exon}s_counts.tsv 


Use RNASEQ-VOILA/scripts/add_phenotype.py with required argument file and optional arguments -t [TYPE] –pval –zscore and -meancorr to set the boundaries for these arguments, default is p-value < 0.01; abs(zscore) > 2.5 and mean corrected > 0. -t 0 is default for outrider input type.

{filename}_pheno_added.tsv is created. This file should be renamed and gzipped to umcu_rnaseq_fib_{untreated/chx}_res_outrider_{gene/exon}s_counts.tsv.gz  
and copied into directory RNASEQ-VOILA/outrider/ 

3.	Add gene associated phenotype to fraser result files:

Use RNASEQ-VOILA/scripts/add_phenotype.py with required argument file and optional arguments -t [TYPE] -pval. Argument -t 2 should be used for fraser files, -pval 0.01 is default.

 
## Manual steps Analyses - Final Master Project Lonneke

For the analysis, featurecounts files are filtered by gene_type=protein_coded. To this extend, an extra column is added to the featureCounts exon_level output. This should be removed for diagnostic runs. To create the reference set, below steps are performed.

1.	In DxNextflowRNA - conf/modules.config - SUBREAD_FEATURECOUNTS_EXON – add gene_type to extraAttributes

```bash
    withName: SUBREAD_FEATURECOUNTS_EXON {
        ext.args = "-g exon_id -M -J --largestOverlap --verbose --extraAttributes gene_id,gene_name,gene_type"
```

2.	To create reference set: use subworkflow featurecounts, for all bamfiles of reference set
3.	copy collected featurecounts files to (reference set) directory and select only gene_type = protein_coding genes/exons 


gene:
```bash
grep -E 'gene_type|protein_coding' ${file} > {refset_directory}/${filename}.{cntrl/chx}.exon.txt
```

exon: 
```bash
grep -E 'gene_type|protein_coding' ${file} | cut -f1-8,10 > {refset_directory}/${filename}.{cntrl/chx}.exon.txt
```
4.	Run OUTRIDER subworkflow: for the analyses I used sample 2021L00022 as ‘sample of interest’ for (cntrl) gene/exon/exon_ratio and (chx) gene versus {reference set MINUS 2021L00022}
5.	HPC: Run outrider separately
```bash
bash run_nextflow_RNAseq.sh {bamfiles directory} {output directory} {email} -entry outrider_entry 
{--inputgene {featurecounts-gene file sample of interest} --refgene {directory featurecounts-gene refset}} 
```

AND / OR

```bash
{--inputexon {featurecounts-exon file sample of interest} --refexon {directory featurecounts-exon refset}}
```

6.	Use EMC APP steps to acquire EMC app files

### EXON RATIO outrider results

1.	Use or create 2 folders for gene and exon (protein coding) counts, only cntrl files (not chx).
2.	Use RNASEQ-VOILA /scripts/exon_ratiofpkm.R to create exon ratio files. Parse arguments for input directories with refsets (featurecounts) 

–refexon: required directory with exon cntrl featurecounts (protein coding)

–refgene: required directory with gene cntrl featurecounts (protein coding)

–output_path: where the exon ratio files are collected
3.	HPC: Create folder for refset and separate (exon ratio) ‘sample of interest’ file 
4.	In NextflowModules/Outrider/1.20.0/outrider.R; function run_outrider,

comment code-block ##FOR EXON LEVEL and comment out code-block ##FOR RATIO FPKM
5.	HPC: Run outrider subworkflow:
```bash
bash run_nextflow_RNAseq.sh {bamfiles directory} {output directory} {email} -entry outrider_entry --inputexon {exon ratio file ‘sample of interest’} --refexon {directory exon ratio refset}
```

6.	metadata and countdata files don’t need to be replaced
7.	EMC APP: In exon ratio output results, replace exon ratio Normalized counts and Mean corrected values with ‘original’ exon-level values from outrider results, with 
RNASEQ-VOILA/scripts/fpkmratio_addexoncounts.R, 
Parse arguments: 

ratioresult (exon ratio  outrider result file)

exonresult (exon level outrider results file)

--output_path
8.	Add phenotypes to Outrider result file: umcu_rnaseq_fib_untreated_res_outrider_exons_counts_countsadded.tsv 

Use RNASEQ-VOILA-MAIN/scripts/add_phenotype.py with required argument file and optional arguments –pval –zscore and -meancorr to set the boundaries for these arguments, default is p-value < 0.01; abs(zscore) > 2.5 and mean corrected > 0.

{filename}_pheno_added.tsv is created. This file should be gzip-ed, renamed to umcu_rnaseq_fib_untreated_res_outrider_exons_counts.tsv.gz  and copied into RNASEQ-VOILA/outrider/ 
