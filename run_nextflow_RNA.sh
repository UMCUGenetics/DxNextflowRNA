#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/projects/RNAseq_Jade/DxNextflowRNA/'

# Set input and output dirs
input_RNA=`realpath -e $1`
bam_boolean=$2
input_WES=`realpath -e $3`
output=`realpath $4`
email=$5
mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name RNAseq
#SBATCH -o log/slurm_nextflow_rna_qc.%j.out
#SBATCH -e log/slurm_nextflow_rna_qc.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

module load Java/1.8.0_60

/hpc/diaggen/software/tools/nextflow run $workflow_path/RNA.nf \
-c $workflow_path/nextflow.config \
--rna_path $input_RNA \
--bam $bam_boolean
--wes_path $input_WES \
--outdir $output \
--email $email \
-profile slurm \
-resume -ansi-log false

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    echo "Remove work directory"
    rm -r work

    echo "Creating md5sum"
    find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

    echo "Data successfully processed to BAM files and VCF files."
    rm workflow.running
    touch workflow.done

    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed
    exit 1
fi
EOT
else
echo "Workflow job not submitted, please check $output for 'workflow.status' files."
fi
