#!/usr/bin/bash

# transcriptome annotation pipeline
# (c) AB 2019
# pure database searches for arbitrary transcriptome files
# Timing is tuned to roughly 30-40k records, naming convention points at transrate "good" contigs
# the argument defines sequence file to be used, check/modify slrun files for details

FILE=$1

echo $FILE
echo 'assuming fasta file is simplified and reasonably non-redundant'
echo 'has .fa extension and is located in the final assembly directory'
#rm $FILE.gff3 $FILE.gff

echo 'submitting DNA-DNA jobs first (infernal nand wise)'
for F in lsu ssu alsu assu blsu bssu tsu MZ12S_p MZ16S_p; do sed s/XXX/$FILE/g slrun_rna_cm.sh | sed s/YYY/$F/g | sbatch; done
for F in ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6; do sed s/XXX/$FILE/g slrun_wise.sh | sed s/YYY/$F/g | sbatch; done 
echo 'submitting single job to generate orfs will wait till it is complete'
# generate orfs
sed s/XXX/$FILE/g slrun_translate_f.sh | sbatch
# wait for it.. Name from slrun_translate_f.sh...
echo 'orf job submitted'
SLOT=0
while [ $SLOT -lt 1 ]; do
 echo 'waiting'
 sleep 3
 SLOT=$(pro-jobs | grep getorf | grep COMPLE| wc -l | awk '{ print $1 }')
done

echo 'submitting blastp -like protein jobs'

for N in {0..3}; do sed s/XXX/$FILE/g slrun_mumsa_f.sh | sed s/NNN/$N/g | sbatch; done
## blast-like jobs... four nodes needed for no more than half an hour
## to search against UniRef90 (2019: about 100mln records)
## because it will take roughly half an hour, should fit plgrid-fast partition
## need further testing for larger transcriptomes (ea. unfiltered data) but seems to work OK

echo 'submit a pyfasta splitting job on translated data for pfam and panther databases'
sed s/XXX/$FILE/g slrun_split.sh | sbatch
# wait for it.. Name from slrun_translate_f.sh...
echo 'splitorf job submitted'
SLOT=0
while [ $SLOT -lt 1 ]; do
 echo 'waiting'
 sleep 3
 SLOT=$(pro-jobs | grep splitorf | grep COMPLE| wc -l | awk '{ print $1 }')
done

# submit python and panther jobs on all chunks
echo 'orf files ready, submitting panther and pfam jobs'
for N in {10..98}; do sed s/XXX/$FILE/g slrun_panther_f.sh | sed s/YYY/$N/g  | sbatch; done
for N in {0..9}; do sed s/XXX/$FILE/g slrun_panther_f.sh | sed s/YYY/'0'$N/g | sbatch; done
for N in {10..98}; do sed s/XXX/$FILE/g slrun_pfam_f.sh| sed s/YYY/$N/g  | sbatch; done
for N in {0..9}; do sed s/XXX/$FILE/g slrun_pfam_f.sh | sed s/YYY/'0'$N/g  | sbatch; done

echo 'all submitted'

