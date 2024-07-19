#!/bin/bash
set -euo pipefail

workflow_path='/hpc/diaggen/software/production/DxNextflowRNA/'

# Set input and output dirs
input=`realpath $1`
output=`realpath $2`
email=$3
optional_params=( "${@:4}" )

mkdir -p $output && cd $output
mkdir -p log

if ! { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
touch workflow.running

sbatch <<EOT
#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name Nextflow_RNASeq
#SBATCH -o log/slurm_nextflow_rnaseq.%j.out
#SBATCH -e log/slurm_nextflow_rnaseq.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --account=diaggen

export NXF_JAVA_HOME='$workflow_path/tools/java/jdk'

${workflow_path}/tools/nextflow/nextflow run \
$workflow_path/main.nf  \
-c $workflow_path/nextflow.config \
--input $input \
--outdir $output \
--email $email \
-profile slurm \
-resume \
-ansi-log false \
${optional_params[@]:-""}
 
if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "RNA workflow completed successfull."
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
