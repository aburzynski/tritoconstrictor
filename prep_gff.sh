#!/usr/bin/bash

# transcriptome annotation pipeline
# (c) AB 2019
# for use after prior DNA (mtDNA and rDNA) database searches.
# Timing is tuned to roughly 30-40k records
# the argument defines sequence file used in searches
# gff with direct DNA hits is generated
FILE=$1

echo $FILE

rm  $FILE.gff

#module add plgrid/tools/python/2.7.14
# make sure the available python can run the script, tested against his version
echo 'collecting DNA (wise) data in gff'

python wise2gff.py $FILE

echo 'augmenting gff with rDNA (infernal)'

# the naming for rRNA profiles may require corrections
for F in lsu ssu alsu assu blsu bssu tsu MZ12S_p MZ16S_p; do python rna2gff.py $FILE'_f.'$F'_cmtbl'; done

echo 'gff ready'

date
