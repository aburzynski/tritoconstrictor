#!/bin/bash -l

#SBATCH -J abus-XXX-YYY
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --time=1:00:00
#SBATCH --output="XXX.YYY_abulog"
#SBATCH --error="XXX.YYY_abuelog"

# example script for submittting read mapping with salmon to a slrum based cluster
# Trinity must be installed
# final assembly must be available under ASSEMBLY, "XXX" will be substituded with assembly name at submission step
# pre-filtered reads for individual samples must be available under READS in paired fq format, "YYY" will be substituted with sample name at submission step.

#export ASSEMBLY=.
#export READS=.
cd $SLURM_SUBMIT_DIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $ASSEMBLY/XXX-trinity/Trinity.fa --seqType fq --left $READS/YYY_1.fq --right $READS/YYY_2.fq --est_method salmon --trinity_mode --output_dir XXX-YYY-abundance-s
