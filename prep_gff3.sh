#!/usr/bin/bash

# transcriptome annotation pipeline
# (c) AB 2019
# for use after prior database searches.
#  the argument defines sequence file used in searches
#  gff3 with protein domains is generated
FILE=$1

echo $FILE

rm $FILE.gff3

#module add plgrid/tools/python/2.7.14
# make sure all dependencies are met in your python environment
echo 'assuming all searches done, just collect the results'
# running chunks in sequence through gff3 generator automagically concatenates results
for N in {0..9}; do python dpfam2gffn.py $FILE.pfam.'0'$N ; done
wc -l $FILE.gff3
for N in {0..9}; do python dpthr2gffn.py $FILE.panther.'0'$N ; done
wc -l $FILE.gff3
for N in {10..98}; do python dpfam2gffn.py $FILE.pfam.$N ; done
wc -l $FILE.gff3
for N in {10..98}; do python dpthr2gffn.py $FILE.panther.$N ; done
wc -l $FILE.gff3
python dblsp2gffn.py $FILE.blastp
wc -l  $FILE.gff3

echo 'gff3 with protein domain matches ready, collect DNA (wise) data in gff separately'
