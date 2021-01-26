#!/bin/bash -l

#SBATCH -J mumsa-NNN-XXX
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 24
#SBATCH --time=01:00:00
#SBATCH --output="XXX.NNN.mumsalog"
#SBATCH --error="XXX.NNN.mumsalog"
#SBATCH -C memfs

module purge

cd $MEMFS

module add plgrid/tools/mmseqs2/11-e1a1c-serial
#export BLASTDB=data/2020
cp $BLASTDB/uniref90_p_NNN_4* $MEMFS
mmseqs createdb $SLURM_SUBMIT_DIR/XXX.orf $MEMFS/XXX.db
mmseqs createindex $MEMFS/XXX.db $MEMFS -s 1 --search-type 2 --compressed 1 
mmseqs search $MEMFS/uniref90_p_NNN_4 $MEMFS/XXX.db $MEMFS/XXX.NNN.out $MEMFS --max-seqs 150 --max-accept 20 -s 1 --threads 24
mmseqs convertalis $MEMFS/uniref90_p_NNN_4 $MEMFS/XXX.db $MEMFS/XXX.NNN.out $MEMFS/XXX.blastp.NNN --format-output "theader,bits,pident,tstart,tend,qheader" --threads 12
cp $MEMFS/XXX.blastp.NNN $SLURM_SUBMIT_DIR
touch $SLURM_SUBMIT_DIR/XXX.NNN.mumsa.OK
# example script for submitting blastp-like search to slurm based cluster
# assumes UniRef90 database split into 4 chunks is present in BLASTDB 
# NNN: database chunk: 0..3
# XXX: proteins in multifasta format with .orf extension, consistent (?) with other db scripts
# creating all intermediate files on MEMFS is desirable, it is a must for IO bound convertalis.
