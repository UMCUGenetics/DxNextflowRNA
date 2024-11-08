[![GitHub Actions CI Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20CI/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/dxnextflowrna/workflows/nf-core%20linting/badge.svg)](https://github.com/UMCUGenetics/dxnextflowrna/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# DxNextflowRNA
UMCU Genetics Nextflow RNA workflow

## Installation
```bash
git clone git@github.com:UMCUGenetics/DxNextflowRNA.git
install.sh
```

## Usage
Note: Replace value of `git_clone_dir` as well as everything between `<>`
```bash
git_clone_dir="/path/to/local/git/clone/"
export NXF_JAVA_HOME="${git_clone_dir}/tools/java/jdk
/hpc/diaggen/projects/woc/rna/tools/nextflow run DxNextflowRNA --input <input_fastq_dir> --outdir <outdir> --email <email>
```

## Citations
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


# Reference files
## GCA_000001405.15 full_plus_hs38d1_analysis_set
The Homo_sapiens GRCh38 GCA_000001405.15 full_plus_hs38d1_analysis_set genome reference is used as main genome reference.
Pre-build star-fusion reference files are used (mainly annotation), which is gencode v22 (GRCh38.p2).
https://github.com/STAR-Fusion/STAR-Fusion/wiki/STAR-Fusion-release-and-CTAT-Genome-Lib-Compatibility-Matrix
https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/

```bash
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play.tar.gz.md5sum
```

## gene bed
Used GTF2BED, Copyright (c) 2011 Erik Aronesty (erik@q32.com)
using perl v5.26.3

## refflat
```bash
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz
```

```bash
singularity shell -B $TMPDIR:$TMPDIR -B /hpc:/hpc https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:357--1
gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons /hpc/diaggen/users/ellen/references/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > /hpc/diaggen/users/ellen/references/GRCh38_gencode_v22_CTAT_lib_Mar012021.ref_annot.gtf.refflat
```

## rRNA intervals
rRNA intervals are retreived from the RSeQC project.
> "We download these ribosome RNAs from UCSC table browser, we provide them here to facilitate users with NO WARRANTY in completeness."

https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_rRNA.bed.gz/download

> RSeQC: quality control of RNA-seq experiments
> _Bioinformatics_ 2012 Aug 15. doi: [10.1093/bioinformatics/bts356](https://dx.doi.org/10.1093/bioinformatics/bts356).

```bash
export NXF_JAVA_HOME="${git_clone_dir}/tools/java/jdk
singularity shell -B $TMPDIR:$TMPDIR -B /hpc:/hpc /hpc/diaggen/software/singularity_cache/broadinstitute-gatk-4.5.0.0.img

gatk CreateSequenceDictionary \
R=/hpc/diaggen/users/ellen/references/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
O=/hpc/diaggen/users/ellen/references/GRCh38_gencode_v22_CTAT_lib_Mar012021.ref_genome.dict
```
- Rename UCSC style to Genbank style; using https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report+ucsc_names.txt
- Exclude all sites that are not included in the references, aka `ALT`, `random` and `chrUn`.

```bash
singularity shell -B $TMPDIR:$TMPDIR -B /hpc:/hpc /hpc/diaggen/software/singularity_cache/broadinstitute-gatk-4.5.0.0.img
gatk BedToIntervalList \
-I /hpc/diaggen/users/ellen/references/hg38_rRNA_genbank_removed_missing_dict.bed \
-O /hpc/diaggen/users/ellen/references/hg38_rRNA_genbank.intervallist \
-SD /hpc/diaggen/users/ellen/references/GRCh38_gencode_v22_CTAT_lib_Mar012021.ref_genome.dict
```

## rRNA database sortmerna
Use rRNA database fasta files provided by sortmerna.
It is recommended to use `smr_v4.3_sensitive_db.fasta` (https://github.com/sortmerna/sortmerna/issues/292)

```bash
wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz

```
File locations are added to `"./assets/sortmerna-db-default.txt"`

Sortmerna is used to remove rRNA from fastq files. Run time can be improved by creating an index first and provide this when running sortmerna a second time.
