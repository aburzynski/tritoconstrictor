#!/bin/bash -l

#SBATCH -J pnthrXXXYYY
#SBATCH --mem-per-cpu=15G
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 5
#SBATCH --time=04:00:00
##SBATCH --output="XXX.YYY.pantherlog"
#SBATCH --error="XXX.YYY.pantherlog"
#SBATCH -C memfs

cd $SLURM_SUBMIT_DIR
cp -L panther.hmm* XXX.orf.YYY $MEMFS
cd $MEMFS
#module add plgrid/tools/hmmer
# make sure the hmmer package is available
hmmsearch --domtblout XXX.panther.YYY --noali --cpu 4 panther.hmm XXX.orf.YYY > /dev/null
cp XXX.panther.YYY $SLURM_SUBMIT_DIR
# example script for submitting panther search to a slurm based cluster
# assumes the presence of links to indexed panther database files in current directory
