#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/users/lonneke/github/DxNextflowRNA'

# Set input and output dirs
input=`realpath $1`
output=`realpath $2`
email=$3

mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

export JAVA_HOME='/hpc/diaggen/software/tools/jdk-18.0.2.1/'  # change java version 
export NXF_JAVA_HOME='/hpc/diaggen/software/tools/jdk-18.0.2.1/'  # change java vesion of nextflow

var_nf=$(sbatch --mem 10G -t 01:00:00 -A diaggen --job-name Nextflow_RNA_testTG -e log/slurm_nextflow_rna_testTG.%j.err \
/hpc/diaggen/software/development/DxNextflowRNA/tools/nextflow run /hpc/diaggen/users/lonneke/github/DxNextflowRNA/test_main.nf  \
-c /hpc/diaggen/users/lonneke/github/DxNextflowRNA/nextflow.config \
--input $input \
--outdir $output \
--email $email)
#-profile slurm \
#-resume -ansi-log false)

echo "${var_nf}"
echo "check directory for output: ${output}"
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi
