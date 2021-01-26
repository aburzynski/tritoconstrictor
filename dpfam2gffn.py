import sys
import re
import os
import os.path

def overlap(x,y):
 # x[0] x.start, x[1] x.end, y[0] y.stat, y[1] y.end
 # single space, no circularity and don't have to care about border cases
 return not (y[1] < x[0] or x[1] < y[0])
 
file_path = sys.argv[1]
# Domtbl - gff3 converter and filter.
# Algo with in-memory dictionary keyed by seqid (complete cds ID) containing a list of lines.
# this may be memory demanding but should lift the major burden from any dowstream parsers
# for panther and pfam gff3 level concatenation is possible since all hits for each orf are in single outfile
# with opening of the name reference file it would be better to run a single instance of the script, not 99!!!
# May run into mamory problem then - copy solution from dblsp2gff script and keep only three best hits
# (on top of score filtering)
# however, for blastp (split database) this is mandatory (same orf can have hits in each outfile)
# This one is pfam-specific.

gff_l={}
with open(file_path, 'r') as f:
 for linia in f:
    if linia[0]<>'#':
     kol = linia.split()
     score = int(float(kol[7]))
     if score >17:
	start = int(kol[19],10)
	end = int(kol[20],10)
        seqid = kol[0]
	gotowa_linia = [start, end, score, ' '.join(kol[22:]).rstrip(),kol[4]] #remember only needed fields

        try:
         lol = gff_l[seqid] #this should rise KeyError if the list is empty
         lsz=len(lol)
         if lsz == 1:
           if lol[0][2] < gotowa_linia[2]:
            gff_l[seqid].insert(0,gotowa_linia)
           else:
            gff_l[seqid].append(gotowa_linia)
         elif lsz == 2:
           if lol[0][2] < gotowa_linia[2]:
            gff_l[seqid].insert(0,gotowa_linia)
           elif lol[1][2] < gotowa_linia[2]:
            gff_l[seqid].insert(1,gotowa_linia)
           else:
            gff_l[seqid].append(gotowa_linia)
         else: 
         #must be 3 then... lets keep it that long... 
         #For panther and unref this is adequate, but for pfam more hits per orf may be good...
         # How to make this part more flexible? Adding fourth element to the list would require
         # 10+ lines of code...
           if lol[0][2] < gotowa_linia[2]:
            del(gff_l[seqid][2])
            gff_l[seqid].insert(0,gotowa_linia)
           elif lol[1][2] < gotowa_linia[2]:
            del(gff_l[seqid][2])
            gff_l[seqid].insert(1,gotowa_linia)
           elif lol[2][2] < gotowa_linia[2]:
            del(gff_l[seqid][2])
            gff_l[seqid].append(gotowa_linia)
        except KeyError:
         gff_l[seqid]=[gotowa_linia]

#print len(gff_l), 'transcripts with pfam hits'

# augmentation of names from "dat" file.
family={}
acc=''
dscrpt=''
with open('pfam_a.dat', 'r') as aux_file:
 contents= aux_file.read()
 records= contents.split('//')
 for r in records:
   lines = r.split('\n')
   del acc, dscrpt
   for ln in lines:
    if ln[:7]=='#=GF AC':
      acc = ln[7:].strip()
    if ln[:7]=='#=GF DE':
      dscrpt = ln[7:].strip()
    if ln[:7]=='#=GF TP':
      tp = ln[7:].strip()
      if tp == 'Domain' :
       if 'protein' not in dscrpt:
        if tp.upper() in dscrpt.upper():
         dscrpt+=' containing protein'
        else:
         dscrpt+=' domain containing protein' 
      elif tp == 'Family':
        if tp.upper() in dscrpt.upper():
         dscrpt+=' member'
        else:
         dscrpt+=' family member'
      break
   try:         
    family[acc]=dscrpt
   except:
    continue    
#print len(family), 'named pfam records in database'

# hierarchical clustering works OK: added best scoring first, then non overlapping worse, 
# limit on acceptable score applied at prefiltering

ftype = 'protein_match'
source = "pfam"
phase = "."
strand ='+'
with open(file_path.split('.')[0] + ".gff3", 'a') as gff:
 for orf in gff_l:
  added=[]
  candidates=sorted(gff_l[orf],key = lambda x: x[2],reverse=True)
  while len(candidates)>0:
   linia = candidates.pop(0)
   if not any(overlap(linia,x) for x in added):
    added.append(linia)
    aux = 'db_xref=PFAM:'+linia[4]+';location='+linia[3]+';description='+family[linia[4]]
    gotowa = [orf,source,ftype,str(linia[0]),str(linia[1]),str(linia[2]),strand,phase,aux]
    gotowa_linia = '\t'.join(gotowa)+'\n'
    gff.write(gotowa_linia)
