#!/bin/bash

#SBATCH --time=24:00:00                  # walltime
#SBATCH --ntasks=1                       # number of tasks
#SBATCH --cpus-per-task=1                # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1                        # number of nodes
#SBATCH -p standard                      # partition(s)
#SBATCH --mem=10G                        # memory per node
#SBATCH -J "Rfam 3D Alignments"          # job name
#SBATCH -o "rfam_3d.out"                 # job output file
#SBATCH -e "rfam_3d.err"                 # job error file
#SBATCH --mail-user=bsweeney@ebi.ac.uk   # email address to message
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -euo pipefail
IFS=$'\n\t'

# This seems to speed up the submission of jobs and thus the overall runtime of
# pipelines.
export NXF_OPTS='-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000'

module load nextflow
nextflow run main.nf -ansi-log false -profile slurm
