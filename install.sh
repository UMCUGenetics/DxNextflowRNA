#!/bin/bash

# Remember the cwd
repo_dir=$PWD

# Create tool directory
mkdir -p ${repo_dir}/tools/

echo "Downloading and installing Nextflow"

export NXF_JAVA_HOME='/hpc/diaggen/software/tools/jdk-18.0.2.1/'
export NXF_VER="25.10.4"


# Download and install nextflow
cd ${repo_dir}/tools/
curl -s https://get.nextflow.io | bash

chmod +x ${repo_dir}/tools/nextflow

${repo_dir}/tools/nextflow -v

echo "Nextflow ${NXF_VER} installed on ${repo_dir}/tools/nextflow"

echo "All done!"
