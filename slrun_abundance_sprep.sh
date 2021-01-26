#!/bin/bash -l

#SBATCH -J prep-XXX-s
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --time=00:15:00
#SBATCH --output="XXX.abuspreplog"
#SBATCH --error="XXX.abusprepelog"

#example script for running abundance analysis, part one: preparing the reference assembly for mapping reads with salmon
cd $SLURM_SUBMIT_DIR
# make sure Trinity is installed
# and the assembly is available
# export ASSEMBLY=.
$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY/XXX-trinity/Trinity.fa --est_method salmon --trinity_mode --prep_reference
