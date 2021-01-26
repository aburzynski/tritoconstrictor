#!/bin/bash -l

#SBATCH -J cmrdna-YYY-XXX
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 10
#SBATCH --time=00:15:00
#SBATCH --output="XXX.YYYlog"
#SBATCH --error="XXX.YYYlog"

cd $SLURM_SUBMIT_DIR

#module add plgrid/tools/infernal
# make sure you have access to the infernal package
# assumes that cm profiles are present in the current directory, "YYY" will be substituted with profile names at submission step
# assumes that the final assembly is available in the ASSEMBLY folder, "XXX" will be substituted with assembly name at submission step
#export ASSEMBLY=.
cmsearch --tblout XXX.YYY_cmtbl --noali --cpu 10 YYY.cm $ASSEMBLY/XXX-trinity/Trinity.fa
