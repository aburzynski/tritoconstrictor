import sys
import re
import os
import os.path
import time

def overlap(x,y):
 # x[0] x.start, x[1] x.end, y[0] y.stat, y[1] y.end
 # single space, no circularity and don't have to care about border cases
 return not (y[1] < x[0] or x[1] < y[0])

def prep_tab_GO(tab_file):
 GO_map={}
 with open(tab_file,'r') as aux_file:
  for ln in aux_file:
   col=ln.split('\t',1)
   key=col[0]
   value=col[1]
   GO_map[key]=value
 return GO_map

def from_tab(tab_file):
 X={}
 print 'reading', tab_file
 with open(tab_file,'r') as aux_file:
  for ln in aux_file:
   col=ln.strip().split('\t')
   key=col[0]
   value=col[1]
   X[key]=value
 return X

big_db={}

dry_run=False
if dry_run:
 mode='w'
else:
 mode='a'
if os.path.isfile('ACC2TaxID.csv'):
 cTax=from_tab('ACC2TaxID.csv')
else:
 cTax={}

t0=time.time()
true_root = sys.argv[1]

## parser of quite specific blastp-like output, defined in mmseq2/convertalis.
## expected order of fields:
## qtitle bitscore pident qstart qend stitle
#   0      1       2      3      4    5
# stitle is from UniRef and must be parsed:
# UnirefclusterID product name, possibly long n=size Tax=taxon name RepID=individualID 
#or more likely
# UnirefclusterID product name, possibly long n=size Tax=taxon name TaxID=number RepID=individualID 
# unfortunately location of the orf in the original fasta is not present in any blastp output
# Therefore it would have to be searched separately from orf file
# However, in mmseq2 the whole defline can be included at will, hence no longer need for this, 
# instead parse the qtitle
# (as in *orf2gff, by passing the relevant fragment to "location") 

# implemented dictionary algo for filter
# all outfiles must be read-in before filtering
# doing it in a loop with guessing of filenames
# must allocate enough memory, but with three elements only kept for each it can be even run directly (mem <5GB)
# further speedup achieved by delaying line processing to the second iteraion
# store raw data - with this memory footprint overhead is of no concern
# have to extract seqid and score and name (name for pre-filtering of GO db)

if len(cTax) == 0:
 file_root=true_root+'.blastpOK.'
else:
 file_root=true_root+'.blastp.'
gff_l ={}
for i in range(4):
 blast_out =file_root+str(i)
 t=time.time()
 print blast_out
 if dry_run:
  d={}
 else:
  try:
   d=prep_tab_GO('uniref2go4db.p_'+str(i))
  except:
   d={}
 print len(d),i,'segment GO terms in', time.time()-t,'s'

 with open(blast_out, 'r') as f:
  for linia in f:
   kol = linia.split('\t')
   score = float(kol[1])
   if score > 17.0:
   #things to store at seqid(0[0]): start(3), end(4), score(1). 
   #Source or precomputed for: pident(2), location(0[1]), name+taxon+description(5). 
        qtitle = kol[0].split(' ',1)
        elements = kol[5].split(' ',1) # at first space to separate the 'name'
        name = elements[0] #name is actually DB reference but with underscore
        name = name.split('_',1)[1]
	kol[5] = elements[1] # store the rest, detailed parsing after passing the filter
        acc = d.get(name,None)
        if acc:
         del d[name] #ensures that assignment is done once only
         big_db[name]=acc
        seqid = qtitle[0]
 	start = float(kol[3])
 	end = float(kol[4])
	gotowa_linia = [start, end, score, kol[2],name,kol[5],kol[3],kol[4],qtitle[1]] 
	#              float  float float pident dbref cluster start end  location-of-orf
	#                0      1    2       3     4     5      6     7        8
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
 t0=t
 t=time.time()
 print t-t0

del d
GO_out=open(true_root+'_uniref2go.filtered', 'w')
print len(gff_l),'orfs with blastp annotations, filtering'
t0=t
ftype = 'protein_match'
source = "uniref90"
phase = "."
strand ='+'

with open(file_root.split('.')[0] + ".gff3", mode) as gff: 
 for orf in gff_l:
  added=[]
  candidates=sorted(gff_l[orf],key = lambda x: x[2],reverse=True)
  while len(candidates)>0:
   linia = candidates.pop(0)
   if not any(overlap(linia,x) for x in added):
        added.append(linia)
        pident = 100*float(linia[3])
        pident=str(pident)
        seqid = orf
	start = linia[6]
	end = linia[7]
        score = str(linia[2])
        stitle = linia[5].rstrip()
        name = linia[4] #name is actually DB reference, pre-extracted
        acc=big_db.get(name,None)
        if acc:
         GO_out.write(name+'\t'+acc)
         del big_db[name]
        aux_d={}
        elements=['',stitle]
        previous_identifier='description=' # This is actually ClusterName
	# complete overhaul, with exhaustive parsing, only RepId and Members are lost (intentional)
	for identifier in ['n=','Tax=','TaxID=','RepID=']:
	 elements=elements[1].split(identifier)
	 aux_d[previous_identifier] = elements[0].strip()
	 previous_identifier = identifier
	aux_d['description=']= aux_d['description='].replace(';',' ').replace('=',' ') 
	# this is needed for poorly formated CusterNames. The presence of  '=' or ';' would ruin parsing of gff3
	aux_d['db_xref=UniProtKB/TrEMBL:']=name
	aux_d['TaxID=']=cTax.get(name,aux_d['TaxID='])
	aux_d['taxon=']=pident[0:4] + '% ' + aux_d['Tax=']
	aux_d['location='] = linia[8]
	del aux_d['Tax=']
	del aux_d['n=']
        aux=';'.join([''.join([k,v]) for k, v in aux_d.items()])
	tbw=[seqid,source,ftype,start,end,score,strand,phase,aux]
	gotowa_linia = '\t'.join(tbw)+'\n'
	gff.write(gotowa_linia)
GO_out.close()
print 'filtered'
