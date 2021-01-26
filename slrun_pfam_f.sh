#!/bin/bash -l

#SBATCH -J pfmXXXYYY
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --time=00:30:00
#SBATCH --error="XXX.YYY.pfamlog"

cd $SLURM_SUBMIT_DIR

#module add plgrid/tools/hmmer
# make sure hmmer is available as well as pfam A hmm file
hmmsearch --domtblout XXX.pfam.YYY --noali --cpu 3 Pfam-A.hmm XXX.orf.YYY > /dev/null

# example script for submitting a pfam search to slurm based cluster
