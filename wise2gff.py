import sys
import re
import os
import os.path
#parser for transcriptome wise files
#iterate over files and append the filtered contents to gff file
#make sure the first column is correct!
# takes one argument: root name for wise/gff files
# it should be possible to extend this to arbitrary wise results, not just mitochondrial

verbose=True
verbose=False
#fasta_file = sys.argv[1]
fasta_name = sys.argv[1]
#mitochondria =['ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5','ND6','ND4L','ATP8']
#mitochondria =['ATP6','CYTB']

mitochondria =['ATP6','COX1','COX2','COX3','CYTB','ND4','ND5']

#linear records assumed with no control on length 

#with open(fasta_file, 'r') as f:
# contents = f.read()
# linie = contents.split('\n')
# seq_length = len(''.join(linie[1:]))

for gene in mitochondria:
 cutoff = 35.0
 if verbose:print 'seeking ', gene
 file_path = fasta_name + '.' + gene+'.gwise'
 while True: #iteration of cutoff if the output is empty not necessarily needed for transcripts but left here...
             #for low scoring output: have to repeat with lower cutoff 
             # if finding score was higher than filtering score due to algorithm pecularities..
  with open(file_path, 'r') as file_input:
	wise_data = file_input.read()
	zakres_linii = wise_data.split("//")
  good = False
  while len(zakres_linii) > 1: #check for empty output
   linie = zakres_linii.pop(0)
   linie = zakres_linii.pop(0)
   if verbose:print linie
   linie = linie.split('\n')
   with open(fasta_name+".gff", 'a') as gff:
    for linia in linie[1:-1]:
        kol = linia.split('\t')
        ftype = kol[2]        
        zapis = ftype == 'cds'
        seqid = kol[0]
        source = kol[1]
	start = int(kol[3])
	end = int(kol[4])
        if not zapis:
          if verbose: print 'supplementary'
          score = kol[5]
        phase = kol[7]
	strand = kol[6]
        if strand == "-":
		start = int(kol[4])
		end = int(kol[3])
        if float(score) < cutoff: #additional filtering, not done by genewisedb, major asset.
          zapis = 0
          if verbose: print 'low score'
	attributes = "product=" + gene 
        gotowa_linia = seqid + '\t' + source + '\t' + ftype.upper() + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes +'\n'
        if zapis:
         if verbose: print 'accepted',seqid,start,end,gene,score,cutoff
         gff.write(gotowa_linia)
         good = True
        else:
         if verbose: print 'rejected', seqid,start,end,gene,score,cutoff
   if verbose:print '.. OK'
   if good:
    if verbose: print 'match found'
#    zakres_linii = []
   else:
    cutoff = cutoff-1
    if verbose: print cutoff, 'was not small enough'
    if cutoff < 10:
      zakres_linii = []
      if verbose: print 'done with this, exiting', gene

  if not good:
   if verbose:print 'not found.'
   break
  else:
    break
