#!/bin/bash -l

#SBATCH -J wiseXXX_YYY
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 12
#SBATCH --time=1:00:00
#SBATCH --output="XXX-YYY.wiselog"
#SBATCH --error="XXX-YYY.wiseelog"

cd $SLURM_SUBMIT_DIR

# make sure your have access to wise2 executables and WISECONFIGDIR points at the correct location
# specify path to hmm profiles for mitochondrial proteins, "YYY" will be substituted with protein name at submission step
# and path to the final assembly folder, "XXX" will be substituted with assembly name at submission step
# the naming convention assumes .fa extension for the final assembly, make sure the file is present
#export MTDNA=.
#export ASSMBLY=.
genewisedb -hmmer -codon codon.table.2 -alg 333 -aalg 333L -pthread -quiet -nohis -sum -gff $MTDNA/YYY.hmm $ASSMBLY/XXX-trinity/Trinity.fa > XXX.YYY.gwise
