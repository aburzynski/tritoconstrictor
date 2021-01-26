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
# for panther and pfam gff3 level concatenation is possible since all hits for each orf are in single outfile.
# This one is panther-specific.

gff_l={}
with open(file_path, 'r') as f:
 for linia in f:
    if linia[0]<>'#':
      kol = linia.split()
      score = int(float(kol[7]))
      if score > 17:
	start = int(kol[19],10)
	end = int(kol[20],10)
        seqid = kol[0]
	gotowa_linia = [start, end, score, kol[3],' '.join(kol[22:]).rstrip()] #remember only needed fields

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
         else: #must be 3 then... lets keep it that long...
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

family={}
with open('names.tab', 'r') as aux_file:
  for ln in aux_file:
    desc = ln.split('\t')
    desc[0]='.'.join(desc[0].split(':'))
    family[desc[0]]=desc[1].rstrip() #should strip newline somehow

#sort by score and save only at most three best hits per orf
ftype = 'protein_match'
source = "panther"
phase = "."
strand ='+'
with open(file_path.split('.')[0] + ".gff3", 'a') as gff: ##!! change w to a in production!!
 for orf in gff_l:
  added=[]
  candidates=sorted(gff_l[orf],key = lambda x: x[2],reverse=True)
  while len(candidates)>0:
   linia = candidates.pop(0)
   if not any(overlap(linia,x) for x in added):
    added.append(linia)
    domain_hit = linia[3].split('.')
    if len(domain_hit) == 3:
          acc = '.'.join(domain_hit[0:2])
    else:
          acc = domain_hit[0]
    ## hack for panther versions mismatch; may be used to limit useless names leaking from panther
    try:
     dscrpt=family[acc]
    except KeyError:
     dom_key=acc.split('.')[0]
     dscrpt=family.get(dom_key,'hypothetical protein') 
    aux = 'db_xref=PANTHER:'+acc+';location='+linia[4]+';description='+dscrpt 
    gotowa_linia = orf+'\t'+source +'\t'+ftype+'\t'+str(linia[0])+'\t'+str(linia[1])+'\t'+str(linia[2])+'\t'+strand+'\t'+phase+'\t'+aux+'\n'
    gff.write(gotowa_linia)
