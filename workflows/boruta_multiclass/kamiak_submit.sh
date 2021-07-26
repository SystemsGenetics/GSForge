#!/bin/bash
#SBATCH --partition=ficklin
#SBATCH --account=ficklin
#SBATCH --job-name=rice
#SBATCH --output=rice%A_%a.out
#SBATCH --error=rice%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0  # Use all the memory of a node.

# Ensure environment variables are configured to use shared resources.
export SINGULARITY_CACHEDIR='/data/ficklin/singularity_cache/'
export NXF_HOME='/data/ficklin/nf_workflows'
export GSFORGE_DEMO_DATA='/data/ficklin/gsforge_demo_data'

# Ensure files created can be accessed by the group.
umask u=rwx,g=rx,o=rx

# Prepare the environment.
module load java nextflow singularity

# Ensure the latest image is being used.
singularity pull -F --name systemsgenetics-gsforge_workflow-latest.img docker://systemsgenetics/gsforge_workflow:latest

# Run the workflow.
nextflow run .main.nf -profile test,kamiak
