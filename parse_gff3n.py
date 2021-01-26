import sys
import re
import os
import os.path
import time
from Bio import SeqIO, SeqUtils, SeqFeature
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import networkx as nx
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
import warnings
warnings.filterwarnings("ignore")
import statistics

overlapl = lambda f1, f2:  any(x in f1 for x in [f2.start, f2.end + (-1)]) or any(x in f2 for x in [f1.start, f1.end + (-1)]) 

rna={'LSUeukarya':'28S large subunit ribosomal RNA','SSUeukarya':'18S small subunit ribosomal RNA','5.8SrRNA':'5.8S ribosomal RNA','12Smitochondria':'12S small subunit ribosomal RNA','16Smitochondria':'16S large subunit ribosomal RNA','LSUbacteria':'23S large subunit ribosomal RNA','SSUbacteria':'16S small subunit ribosomal RNA','LSUarchaea':'23S large subunit ribosomal RNA','SSUarchaea':'16S small subunit ribosomal RNA'}

k2i={'bacteria':'2','archea':'2157','archaea':'2157','eukarya':'2759','mitochondria':'2759','rRNA':'2759'}    
def diminute_end(loc):
 x = loc
 xe = x.end - 1
 x = FeatureLocation(x.start, xe, x.strand)
 return x

def diminute_start(loc):
 x = loc
 xs = x.start + 1
 x = FeatureLocation(xs, x.end, x.strand)
 return x

def loc2gb(loc):
    from StringIO import StringIO
    out_handle = StringIO()
    f = SeqFeature.SeqFeature(loc,'misc',{})
    rec = SeqRecord(r.seq)
    rec.features.append(f)
    out_handle.write(rec.format('gb'))
    result = out_handle.getvalue().split('misc')[1].split('ORIGIN')[0].strip()
    return result

def gb2loc(location_text): #expected input format: pos:gb_location
 a=location_text.split(':')
 s=1
 loc = FeatureLocation(0,sequence_len, strand=s)
 if a[0] == 'pos':
    if a[1][:10] == 'complement' :
      s = -1
      a[1] = a[1][11:-1]
    p = a[1].split('..')
    if len(p) == 1:
      p.append(p[0])
    loc = FeatureLocation(int(p[0])-1,int(p[1]),strand=s)
 return loc

def make_gff_list(gff_name): # actual gff3 parser - for protein hits
 disqualifying= ['NOT NAMED','UNCHARACTERIZED','GENOME SHOTGUN','UNNAMED','HYPOTHETICAL','EMB|','HIGHLY SIMILAR TO DELETED','GENE RICH CLUSTER','GENE REGULATED BY']
 undesired=['PLASMODIUM', 'TRUNCATED','PSEUDOGENE','FRAGMENT','LOW QUALITY PROTEIN','HOMO SAPIENS','HUMAN','DROSOPHILA','YEAST','S. CEREVISIAE','SACCHAROMYCES CEREVISIAE','S. CEREVISIASE','XENOPUS','S. POMBE','SCHIZOSACCHAROMYCES POMBE','E. COLI', 'ECOLI', 'ARABIDOPSIS','CHLAMYDOMONAS','H. INFLUENZAE','H. SAPIENS','MOUSE','C. ELEGANS','PROTEIN PROTEIN','CONSERVED PROTEIN', 'PROBABLE ','PREDICTED PROTEIN','CONSERVED ','EXPRESSED PROTEIN','OPEN READING FRAME',' GENE','GENE/','COPIES','GENOME DUPLICATE','ISOFORMS ','PARTIAL']
 acceptable={'Homo B':'Homo-B','granulins':'granulin','Granulins':'Granulin','transcriptases':'transcriptase','neuropeptides':'neuropeptide','Staphylococcal nuclease':'SNase','Staphylococcal':'STAPH','rat sarcoma':'sarcoma','domains':'domain','pathways':'pathway','motifs':'motif','repeats':'repeat','hands':'hand','cells':'cell','pigments':'pigment','muscles':'muscle','.':' ','ortholog':'-like','paralog':'-like','homologous':'similar to','homolog':'-like','homologue':'-like','orthologue':'-like','orphan protein':'protein','novel':'putative','Novel':'Putative',':':' ','_':' ','(':' ',')':' ','[':' ',']':' '}
 dspace={'  ':' '}
 gff_dict ={}
 with open(gff_name, 'r') as gff_input:
  for raw_gff in gff_input:
   gff = raw_gff.split('\t')
   if len(gff) > 8:
    start = int(gff[3])
    end = int(gff[4])
    f_type = gff[2]
    qual = {}
    key_val_list = gff[8].split(';')
    for key_val in key_val_list:
     key_list = key_val.split('=',1)
     if len(key_list)>1:
      qual.update({key_list[0]:key_list[1].rstrip()})
    direct = None
    loc_str = qual['location']
    loc_part = loc_str.split()
    s_nt = int(loc_part[0][1:])-1
    if len(loc_part) <4: 
      direct = 1
      x_nt = (start-1)*3 + s_nt
      y_nt = (end)*3 + s_nt
    else:
      direct = -1
      x_nt = s_nt - 3*(end) +1
      y_nt = s_nt - 3*(start-1) +1
    qual.update({'source':gff[1], 'score':float(gff[5])})

    # additional filtering of "description"
    # remove brackets and certain undesired phrases (fragment, low quality protein etc.)
    # consider either more or less sophistication as this gets you nowhere; 
    # still 1000s of products don't meet ncbi guidelines

    product = qual.get('description','')
    if product:
     #remove entire description
     if any(x in product.upper() for x in disqualifying):
        del qual['description']
#        print product,'R'
     else:
      #remove phrases
      product_up=product.upper()
      fixes=0
      for x in undesired:
       while x in product_up:
         fixes += 1
         idx = product_up.find(x)
         product = product[:idx]+product[idx+len(x):]
         product_up = product.upper()
      #cut and remove bogus isoform invocation
      if 'isoform' in product:
         idx=product.index('isoform')
         product=product[:idx]
         fixes += 1
      #substitute some other phrases with more "politicaly correct" versions
      for x in acceptable:
       while x in product:
        fixes += 1
        idx = product.find(x)
        product = product[:idx]+acceptable[x]+product[idx+len(x):]
      #cut at "-like" (keed "like..." though - investigate)
      if '-like' in product:
         idx=product.index('like')
         product=product[:idx+4]
         fixes += 1
      #final cleanup of doublespaces
      for x in dspace:
       while x in product:
        fixes += 1
        idx = product.find(x)
        product = product[:idx]+dspace[x]+product[idx+len(x):]
      #tbl2asn silently substitutes "illegal" trailing characters with underscores and then complains! Beware!
      qual['description'] = product.strip().strip(",-/' ")
      if qual['description'].upper() in ['PROTEIN','RELATED',''] or fixes > 4:
        del qual['description']

    if x_nt <0:
     x_nt = BeforePosition(0)
    else:
     x_nt = ExactPosition(x_nt)
    y_nt=ExactPosition(y_nt)
    tmp = gff[0].split('_')
    transcript = '_'.join(tmp[:-1])
    qual['cds_id'] = tmp[-1]
    loc = FeatureLocation(x_nt,y_nt,strand=direct)
    clean_feature = SeqFeature.SeqFeature(loc, f_type, {})
    clean_feature.qualifiers = qual
    try:
      gff_dict[transcript].append(clean_feature)
    except KeyError:
      gff_dict[transcript]=[clean_feature]
#   print '.'
 return gff_dict # now the dictionary stores lists of annotations ready to be combined with records

# for direct, nt level annotations (ea. rna)
def make_gff_list_rna(gff_name):
 gff_dict={}
 with open(gff_name, 'r') as gff_input:
  for raw_gff in gff_input:
   gff = raw_gff.split('\t')
   if len(gff) > 8:
     start = int(gff[3])
     end = int(gff[4])
     f_type = gff[2]
     direct = None
     if gff[6] == '+': 
      direct = 1
     elif gff[6] == '-':
      direct = -1
     qual = {}
     key_val_list = gff[8].split(';')
     for key_val in key_val_list:
      key_list = key_val.split('=',1)
      qual.update({key_list[0]:key_list[1].rstrip()})
     qual.update({'source':gff[1], 'score':float(gff[5])}) 
     location = FeatureLocation(start-1, end, strand=direct)
     new_feature = SeqFeature.SeqFeature(location, type=f_type, qualifiers=qual)
     transcript = gff[0]
     try:
      gff_dict[transcript].append(new_feature)
     except KeyError:
      gff_dict[transcript]=[new_feature]
 return gff_dict

#overlap-based clustering. New, safer but a bit clumsy... Check gff2gb...
def rclust(L1):
 LL =[]
 while len(L1) > 0:
  x = L1.pop(0)
  candidates = [x]
  while True:
   match = []
   for f in candidates:
    for g in L1:
     if overlapl(f.location,g.location):  
      match.append(g)
   for g in match:
    if g in L1:
     L1.remove(g)
    if g not in candidates:
     candidates.append(g)
   if match == []:
    break
  LL.append(candidates)
 return LL #list-of-lists 

# CDS cluster construction
# pop sequentially from domains and check all cds_ids
# construct list-of-lists for each cds_id
# make CDS feature for each cluster

def cluster(cds): #cluster by cds_id which is based on getorf "names" stored in gff
 d={} #d is dictionary of feature lists keyed by cds_id
 while len(cds) >0:
  f=cds.pop()
  p = f.qualifiers['cds_id']
  del f.qualifiers['cds_id']
  try:
    d[p].append(f) #creating list-of-features for each cds separately
  except KeyError:
    d[p]=[f]
 #iterate over cds_ids, create CDS feature for each and store it together with the respective list in d
 for cds_id in d:
  dom=sorted(d[cds_id], key = lambda x: x.qualifiers['score'], reverse=True)
  loc_str = dom[0].qualifiers['location'] #all elements have the same location (unavoidable redundancy)
  for f in d[cds_id]:
   del f.qualifiers['location'] # this operates on dom elements actually (alias)
  loc_part = loc_str.split()
  if len(loc_part) <4: 
      direct = 1
      x_nt = int(loc_part[0][1:])-1
      y_nt = int(loc_part[2][:-1])
  else:
      direct = -1
      y_nt = int(loc_part[0][1:])
      x_nt = int(loc_part[2][:-1])-1
  loc = FeatureLocation(x_nt,y_nt,strand=direct)
  cds_f = SeqFeature.SeqFeature(loc, 'CDS', {})
  cds_f.qualifiers['score'] = dom[0].qualifiers['score'] # best score... but not necessary the best idea...
  references=[]  
  for each in dom:
   if float(each.qualifiers['score']) > 20: #consider tuning this up... to avoid weak hits leaking into submission
    references.append(each.qualifiers['db_xref'])
  if len(references)>0:
   tset=set()
   fset=set()
   empty=set()
   key_set=[]
   for each in references:
    component=each.split(':')
    database = component[0]
    key=component[1].split('.')[0]
    if database in ['PANTHER','PFAM']:
     key_set.append(key)
     if database == 'PANTHER':
      pthr_points_at=pth.get(key,empty)
      tset.update(pthr_points_at)
     else:
      pfam_points_at=pfm.get(key,empty)
      fset.update(pfam_points_at)
   key_set=set(key_set)
   iprset = tset | fset
   if len (iprset) > 1: #if more possibilities do intersection, should fish out additional matches
    iprset= tset & fset
   if len(iprset) == 1:
     #profile points at one ipr entry only, common (intersection) or individual (union) 
     # cases like common pfam domains pointing at many iprs are not considered
     # check if set of acc nos for this particular ipr matches the annotation (set of database keys)
     ipr_candidate=iprset.pop()
     dom_set=ipr_d[ipr_candidate] #if we are here then it is guaranteed to exist, no need to secure the call with get
     if key_set == dom_set:
      references.append('InterPro:' + ipr_candidate)    
   cds_f.qualifiers['db_xref']=list(set(references))
  dom.append(cds_f)
  d[cds_id]=dom #after this first has the best scoring domain, last is cds.
 return d

def parse_tbl(tbl_f):
 d={}
 with open(tbl_f,'r') as fh:
  #read the whole contents
  tbl = fh.read()
 entrees=tbl.split('>Feature')
 for datalines in entrees:
  line=datalines.split('\n')
  offset = line[0].find('TRINITY')
  transcript = line.pop(0)[offset:].rstrip()
#  print transcript
  while len(line) > 0:
   flocation = line.pop(0)
   element=flocation.split('\t')
   loc_str= '-'.join(element[:2])
   while len(line) > 0:
     fqualorig = line.pop(0)
     fqual = fqualorig.split('\t')
     if fqual[0] <> '':
      line = [fqualorig] + line #run into another location, re-create and break the inner loop.
      break
     elif len(fqual) == 1: # last line! splitting worked strangely...
      
      continue
     if fqual[3]=='product' and element[2]=='CDS':
      d[transcript+':'+loc_str]=fqual[4].strip()
#      print transcript+':'+loc_str+'\t'+d[transcript+':'+loc_str]
 return d

def prep_GO(pthr_go): # return one dictionary of flat lists (of GO term numbers), keyed by panther ID
 GO_map={}
 with open(pthr_go, 'r') as aux_file:
  for ln in aux_file:
    desc = ln.split('\t')
    if len(desc) > 2:
     key_p=desc[0].replace(':','.')
     go_list=[]
     go_part=desc[2].split(';')
     for each in go_part:
      number=each.split('#GO:')
      if len(number) >1:
       go_list.append(number[1])
     go_part=desc[3].split(';')
     for each in go_part:
      number=each.split('#GO:')
      if len(number) >1:
       go_list.append(number[1])
     go_part=desc[4].split(';')
     for each in go_part:
      number=each.split('#GO:')
      if len(number) >1:
       go_list.append(number[1])
     if len(go_list)>0:
      GO_map[key_p]=go_list
 return GO_map
 
def read_GO(file_root): # blindly read-in table of GO sets into dictionary
 fname=file_root+'.filtered'
 d={}
 with open(fname,'r') as fh:
  for ln in fh:
   col=ln.strip().split('\t')
   key=col[0]
   value=col[1].split(',')
   d[key]=set(value)
 return d

def from_tab(tab_file): # general procedure for reading dictionary from tab delimited tsv
 d={}
 with open(tab_file,'r') as fh:
  for ln in fh:
   col=ln.strip().split('\t')
   key=col[0]
   value=col[1]
   d[key]=set(value)
 return d
  
def prep_IPR_GO(ipr_go): # return similar structure as the above but keyed by ipr acc no
 GO_map={}
 with open(ipr_go, 'r') as aux_file:
  contents=aux_file.read()
 lines=contents.split('\n')
 for ln in lines[3:]:
  col=ln.split(':')
  if len(col) > 1:
   ipr_key=col[1][:9] # interpro accessions are always 9 digit long
   go_term=col[-1].strip()
   try:
    GO_map[ipr_key].append(go_term)
   except KeyError:
    GO_map[ipr_key]=[go_term]
 return GO_map
 
def classify(obo_file): 
# return a DiGraph with GO classification, along with dictionaries of names and alternative identifiers
# version for filtered databases, should not have non-is_a conncetions
# but database for filtering should have them
# store these as an additional set of edges
# store also filter of terms not allowed in annotations
# these are needed (?) because filtering at consolidation stage is non-trivial.
 singular_terms=['id','name','namespace','def','is_obsolete']
 N_d={}
 A_d={}
 BAD=[]
 fedge_list=[]
 with open(obo_file,'r') as aux_file:
  contents=aux_file.read()
 terms=contents.split('\n\n')
 node_list=[]
 edge_list=[]
 for candidate in terms:
  line=candidate.split('\n')
  if line[0].strip()=='[Term]':
   d={}
   for each in line[1:]:
     col=each.split(':',1)
     try:
      d[col[0]].append(col[1])
     except KeyError:
      d[col[0]]=[col[1]]
   k_list = list(d)
   for x in k_list:
    if x in singular_terms:
     d[x]=d[x][0]
   is_obsolete=d.get('is_obsolete','').strip()
   if not (is_obsolete == 'true'):
    term=d['id'].split(':')[1].strip()
    is_a=d.get('is_a',[])
    N_d[term]=d['name']
    alternative_ids=d.get('alt_id',[])
    for x in alternative_ids:
     alt=x.split(':')[1].strip()
     A_d[alt]=term
    subsets=d.get('subset',[])
    for x in subsets:
     if 'gocheck' in x:
      BAD.append(term)
    
    rlshp=d.get('relationship',[])
    for relation in rlshp:
     if 'part_of' in relation:
      term2=relation.split(':')[1][:7]
      edge=(term,term2)
      fedge_list.append(edge)
    
    node_list.append(term)
    for relation in is_a:
     term2=relation.split(':')[1][:7] # assuming GO identifiers are 7 char long, better use more intelligent splitting
     edge=(term,term2)
     edge_list.append(edge)
 GO=nx.DiGraph()
 GO.add_nodes_from(node_list)
 GO.add_edges_from(edge_list)
 return GO,N_d,A_d,BAD,fedge_list

 # reduce and filter the mapping to
 # non-redundant non-obsolete and consistent (non-alternative)
 # terms, stored as sets in a dictionary keyed by original keys (accessions)
 ##Assumes global variables:##
 # alt: dictionary of alternative ids 
 # GO: graph of GO term 'is_a' relations
def filter_map(original_map):
 filtered_map={}
 for acc in original_map:
  all_desc=set()
  for term in original_map[acc]:
   f_term=alt.get(term,term)
   try:
    tset=nx.descendants(GO,f_term)
   except:
    tset=set([f_term]) #if no descentants, remove the term itself - either it is obsolete or trivial
   all_desc.update(tset)
  new_list = [] 
  for term in original_map[acc]:
   f_term=alt.get(term,term)
   if f_term not in all_desc:
    new_list.append(f_term)
  if len(new_list)>0:    
   filtered_map[acc]=set(new_list)
 return filtered_map
 
 # similar procedure for a single group of arbitrary GO terms 
 # x can be a set or list of GO terms. GO is taken from global. 
 # only for confirmed terms, no error handling!
 # returns list... If x was a set then this should be OK, otherwise make sure its non-redundant.
def filter_list(x):
 y=[]
 bad=set(BAD)
 for term in x:
  tset=nx.descendants(GO,term)
  bad.update(tset)
 for term in x:
  if term not in bad:
   y.append(term)
 return y

## For taxonomy based filtering
def get_tax(tab_file): 
 d={}
 edge_list=[]
 node_list=['1']#'root' node added manually
 with open(tab_file,'r') as aux_file:
  for ln in aux_file:
   col=ln.strip().split('\t')
   node_list.append(col[0])
   edge_list.append((col[0],col[1]))
   d[col[0]]=col[2]
 TAX=nx.DiGraph()
 TAX.add_nodes_from(node_list)
 TAX.add_edges_from(edge_list)
 return TAX,d

def revname(d): # reverse dictionary. could reverse any dictionary with hashable values. Produces lists of keys.
 n_d={}
 for tax_id in d:
  try:
   n_d[d[tax_id]].append(tax_id)
  except KeyError:
   n_d[d[tax_id]]=[tax_id]
 return n_d

def rca(tid_list):
  cu_set = set(tid_list)
#  tid = str(sorted([int(x) for x in tid_list])[0])
  tid = tid_list.pop()
  a_set = nx.ancestors(TAX,tid)
  a_set.add(tid)
  while not cu_set <= a_set:
   p_l = list(TAX[tid])
   tid = p_l.pop()
   a_set = nx.ancestors(TAX,tid)
   a_set.add(tid)
  return tid

###START HERE

t1=time.time()
fasta_name = sys.argv[1]
print 'reading', fasta_name
del_dict={}
gff_name = sys.argv[2]
name = gff_name.split(".")
output_name = name[0] + ".gb"

# read all records into memory
records = list(SeqIO.parse(fasta_name, "fasta"))
t2=time.time()
print t2-t1,'s\n', len(records), 'records in, parsing gff3'

# read all annotations also into memory. For typical size transcriptomes this is not a problem.
gff_l = make_gff_list(gff_name)

rdict = SeqIO.to_dict(records)
for f in gff_l:
 rdict[f].features=gff_l[f]

t3=time.time()
print t3-t2, 's\n', len(gff_l),'records with protein hits, preparing aux data'

pth={}
pfm={}
ipr_d={}
ipr_file='interpro.tab'
with open(ipr_file,'r') as ipr:
 for ln in ipr:
  col=ln.strip().split('\t')
  key = col[2].split('.')[0]
  value = col[0]
  if col[1] in ['PANTHER','PFAM']:
   try:
    ipr_d[value].append(key)
   except KeyError:
    ipr_d[value]=[key]
   if col[1]=='PANTHER':
    try:
     pth[key].append(value)
    except KeyError:
     pth[key]=[value]
   else:
    try:
     pfm[key].append(value)
    except KeyError:
     pfm[key]=[value]   
for each in pth:
 pth[each]=set(pth[each])
for each in pfm:
 pfm[each]=set(pfm[each])
for each in ipr_d:
 ipr_d[each]=set(ipr_d[each])

print ' pfam accessions with interpro reference:',len(pfm)
print ' panther accessions with interpro reference:',len(pth)

dup_dict={}
filtered_fasta=name[0]+'.filtered'
# this should be done on sets, not dictionaries
if os.path.isfile(filtered_fasta):
 print ' detected filtered file, processing.'
 print ' Consider re-running dbsearch with filtered file instead.'
 good_records=list(SeqIO.parse(filtered_fasta,"fasta"))
 gdict = SeqIO.to_dict(good_records)
 all_set = set(list(rdict))
 good_set = set(list(gdict))
 dup_set = all_set - good_set
 for r in dup_set:
    dup_dict[r]=True
    rdict[r].description = 'duplicate'
 del gdict, all_set, dup_set, good_records
 print len(dup_dict), ' records marked as duplicates'

print ' Preparing GO terms.'
rGO,ref,alt,BAD,fedges = classify('go-basic.obo')
fedges.append(('0005622','0005623')) # missing connection between 'cell' and 'intracellular', needed for proper filtering
fadd=nx.DiGraph()
fadd.add_edges_from(fedges)
GO=nx.compose_all([rGO,fadd])
del fedges
 ## Prepare the GO classification sets here.
 ## This is for (presumabely) four sets of terms
 ## zero - master, the root set of top terms (presumabely just three of them)
 ## 1 - from root 0, 2 -from root 1, 3 from root 2 etc.
 ## these level-less sets under each root term, consist of all "ancestors" of each root term
 ## "level" parameter could be calculated later as the distance from/to the root term, as needed.

ds=set()
for term_g in nx.attracting_components(rGO):
  term=term_g.pop()
  ds.add(term)
classify_list=[ds]
class_name=['root']

for root in classify_list[0]:
  ds=nx.ancestors(rGO,root)
  classify_list.append(ds)
  class_name.append(root)
classify_d=dict(zip(class_name,classify_list))

print ' understand',len(GO),'GO terms in categories:'
for ctgry in classify_d['root']:
 print '\t',ref[ctgry], len(classify_d[ctgry])

if os.path.isfile('interpro2go.filtered'):
 GO_ipr=read_GO('interpro2go')
else:
 ipr_map=prep_IPR_GO('interpro2go')
 GO_ipr=filter_map(ipr_map)
print ' interpro accessions with GO terms:',len(GO_ipr)
## InterPro accessions are the cleanest but they cover usually (?) ca. 10% of transcripts only

if os.path.isfile('panther2go.filtered'):
 GO_pth=read_GO('panther2go')
else:
 pthr_map=prep_GO('PANTHER15.0_HMM_classifications')
 GO_pth=filter_map(pthr_map)
print ' panther accessions with GO terms:', len(GO_pth)
### PANTHER GO assignments tend to be inconsistent
### there is no way to clean this up automatically and fully reliably
### Leave the code for now but use with caution

uniroot=name[0]+'_uniref2go'

if os.path.isfile(uniroot+'.filtered'):
  GO_uni=read_GO(uniroot)
else:
  try:
   GO_uni=read_GO('uniref2go')
  except:
   GO_uni={}
print ' UniRef90 accessions with GO terms:', len(GO_uni)
## UniRef GO dictionary is HUGE. Would need +20GB RAM to use the whole db. 
## Therefore only using of prefiltered subsets is implemented for this database.
## Either the whole UniRef90 (still large) or transcriptome-specific
## table filtered with blast results (should be managable)
## generation of the filtered db at dblsp2gff stage is OK, overhead do not exceed 30s
## However, the split and prefiltered files must be prepared in advance
## They are not strictly version-specific, 
## the newest assignments can be used on old databases (but the reverse is risky)

itaxon={} 
#keys will be taxons encountered, values will be lists of records. 
#TAX,i2n=get_tax('taxonomy.tsv')
TAX,i2n=get_tax('u90-simplified.tsv') #version adequate for uniref90 flavor
#TAX,i2n=get_tax('simplified.tsv')
#n2i=revname(i2n)
#actually, n2i may not be needed anymore... Check.

#pre-calculate "ancestors" of predefined "kingdoms" for rrna-based taxonomy 
# these are de facto lists of all taxons contained in a kingdom
# pre-calculating these takes memory but the speedup is needed
k2anc={}
danc={}
idset=set()
for k in k2i:
 idset.add(k2i[k])
for tid in idset:
 danc[tid] = nx.ancestors(TAX,tid)
 danc[tid].add(tid)
for k in k2i:
 k2anc[k]=danc[k2i[k]]
del danc, idset 

print ' taxons from NCBI database:', len(TAX)
# TAX database: DiGraph similar to GO, with one root, contains TaxIDs 

t4=time.time()
print t4-t3, 's\naggregating annotations'

len_dict={}

for fr in rdict:
    rec=rdict[fr]
    len_dict[fr]=len(rec)
    rec.seq.alphabet = IUPAC.ambiguous_dna
    ttable=CodonTable.unambiguous_dna_by_id[1]
    stop_codons=ttable.stop_codons
    start_codons=ttable.start_codons
    rec.name = rec.name[10:] #don't need that for imgt format but for gb it is a dirty workaround for LOCUS parsing problem
    if len(rec.features) > 0:
      clust = cluster(rec.features) #note that this removes all features from rec as a side effect.
      candidates=[]
      candidate_cds=[]
      for cds_id in clust:
        f_set=clust[cds_id]
        clean_feature = f_set[-1]
        f_set = f_set[:-1]
        #transfer the best scoring description as product to the CDS, but only if meaningful
        ## this is the most problematic part of the entire pipeline##
        for f in f_set:
         try:
          clean_feature.qualifiers['product'] = f.qualifiers['description'].strip()
#          if f.qualifiers['source'] in ['panther','pfam']:
          break
         except:
          clean_feature.qualifiers['product'] = 'hypothetical protein'
          continue
        tid_list=[]
        for f in f_set: 
          organism = f.qualifiers.get('taxon','0%').split(' ',1)
          identity = float(organism[0].strip('%'))
          score = f.qualifiers.get('score',0)
          if identity > 85 and score > 40: 
          # bogus taxons give a bad clue only so just ignore them
          # however, for smaller samples lowering this may be beneficial as whitelisting would be based on stronger input
          # maybe differentiate whitelisting and filtering steps? But how?
           tid = f.qualifiers['TaxID']
           tid_list.append(tid)
        possibilities=len(tid_list)   
        if possibilities>1:
         keep=[]
         for tid in tid_list:
          if tid in TAX:
            keep.append(tid)
         if len(keep)>1:
          tid=rca(keep)
         elif len(keep) == 1:
          tid=keep.pop()
         else:
          possiblities=0             
        if possibilities > 0 and tid in TAX:
             orgnsm=i2n[tid]
             clean_feature.qualifiers['taxon'] = orgnsm
             clean_feature.qualifiers['TaxID'] = tid

        loc=clean_feature.location
# dont't have to fix start codons because getorf would not generate such cases. Perhaps this is actually bad?
# NEED to fix all orfs to include stop codon!
        if clean_feature.strand == -1:
             new_start=loc.start-3
             if new_start < 0:
               clean_feature.location = FeatureLocation(BeforePosition(0),loc.end,strand=-1)
             else:
               clean_feature.location = FeatureLocation(new_start,loc.end,strand=-1)
        else:
             new_end=loc.end+3
             if new_end > len(rec.seq):
               clean_feature.location = FeatureLocation(loc.start,AfterPosition(len(rec.seq)),strand=1)
             else:
               clean_feature.location = FeatureLocation(loc.start,new_end,strand=1)

        # quality filtering? Overlap check? Anything to avoid bogus annoatations? 
        # need to store the whole structure temporarily (new list), 
        # then could sort/filter the list and add items to rec.features.
        # cds and good orf counting could be done in the same loop
        candidate_cds.append(clean_feature)
        f_set.append(clean_feature) #cds is again the last one
        candidates.append(f_set) # list-of-lists
        
      cds_set=rclust(candidate_cds) 
      # cluster by overlap. Each cds_set is a list of cds features.
      # not sure how the list is ordered, but it was sorted at some point...
      # just to be on the safe side sort it again and prefilter
      good_cds=[]
      for each in cds_set:
       best=sorted(each, key = lambda x: x.qualifiers['score'], reverse=True)[0]
       if best.qualifiers['score'] > 17:
        good_cds.append(best)

      accepted = 0
      good_orf_count = 0    
      longest_orf = 0
      if len(good_cds) > 0:
       candidates = sorted(candidates, key = lambda x: x[-1].qualifiers['score'], reverse=True)
       for f_set in candidates:
         cds = f_set[-1]
         if cds in good_cds:     
             long_orf=len(cds)
             if long_orf > longest_orf:
               longest_orf = long_orf
             if cds.qualifiers['score'] > 20: 
               rec.features.extend(f_set)
               accepted+=1
               good_orf_count+=1
             elif accepted <2: 
             # after prefiltering this would never be used... so prefilter is more relaxed.
               rec.features.extend(f_set)
               accepted+=1
# store accepted and good_orf_counts in each record in dbxrefs list...
# lets try... seems to work but only internally. 
# Really dirty hack, consider using some custom field instead (it is saved but not imported anyway)
# or indeed a custom dictionary with the same key could be used to store these internally and save to csv.
# But saving somewhere in the record structure would be preferable (could use the same for replacement records)
# -for now this does not work for them-
      rec.dbxrefs=[str(good_orf_count),str(longest_orf)] 
      try:    
       longest_good_feature = sorted(good_cds, key = lambda x: len(x)*x.qualifiers['score'], reverse=True)[0]
       descript = longest_good_feature.qualifiers['product']
       if rec.description == 'duplicate':
         descript += ' duplicate'
       rec.description = descript
      except:
       pass
      ## Sophisticated taxon assignment IMPLEMENTED! Two things: 
      ## -uses TAX database (and taxids) - would require modification of gff3 generator - DONE
      ## -make sure the best possible id is passed, must come from the best scorer, 
      ##  must be high% but dosen't have to come form the best feature (?)
      ##  this accounting/aggregation should be done much earlier at CDS definition? - done, OK
      ##  then here just would have to check if best CDS has taxon designation and pass it to record
      ##  as well as store in db.
      ##  MOVE to taxid based classification completed, and uniref TaxID bug fixed.
      ##  Reliable TaxID stored in itaxon database. Names stored in record, but this is never actually used for filtering.
      ##  ONE more caveat: the 80% limit may be actually too relaxed for some very conservative
      ##  genes. Any firm taxonomy based on such genes would require nearly perfect match and high score.
      ##  Consider making a list of candidates either entirely excluded or requiring better support.
      ##  Such a list could be prepared after inspecting deletion files from several assemblies.
      
      try:
          orgnsm=longest_good_feature.qualifiers['taxon']
          rec.annotations['organism'] = orgnsm
          tid=longest_good_feature.qualifiers['TaxID']
          rec.annotations['TaxID']=tid
          # store list of identifiers. Avoid counting the same multiple times (sets)
          try:
           itaxon[tid].add(fr)
          except KeyError:
           itaxon[tid]=set([fr])
      except:
        pass
      # all features added, transcript annotated. Augment it with GO terms here.
      # collect possible lists based on db_xref
      # d.get() should work nicely around missing terms
      # there are currently 3 dictionaries: GO_pthr and GO_ipr and GO_uni
      # all non-redundant mappings (set) end up in one place 
      # store it with a cds feature as a custom qualifier
      
      for each in rec.features:
       if each.type=='CDS':
         go_terms=set()
         try:
          dbreflist=each.qualifiers['db_xref']
         except:
          dbreflist=[] 
         for element in dbreflist:
           db_key=element.split(':')[1]
           go_terms.update(GO_ipr.get(db_key,[]))
           go_terms.update(GO_pth.get(db_key,[]))
           go_terms.update(GO_uni.get(db_key,[]))
         if len(go_terms) >0:  
           each.qualifiers['GO_terms']=filter_list(go_terms)
#removing redundancy should be enough if annotations are consistent and databases filtered
# BUT it should be better to filter using the "aggressive" database with part_of connections and gocheck filters
# IF descendants filtering should be desirable
# YES, seems to be needed. Check the results!!!

      if good_orf_count == 1:
        rec.description +=', simple transcript'
      else:
        rec.description +=', problematic transcript, '+str(good_orf_count)+' orfs'
      
print 'aggregated',
t6=time.time()
print t6-t4,'s\nparsing additional dna-level gff'


# integrated DNA annotations
# read also all DNA annotations into memory; there usually are very few of them...
# currently these are mainly used for marking records for deletion

gff_r = make_gff_list_rna(name[0]+'.gff')

id_file= name[0]+'.id'
eliminate_mt = os.path.isfile(id_file)
del_mt={}
if eliminate_mt:
 print 'identified mitotranscripts will be removed, may add them back by editing id file'

rnid_file=name[0]+'.rnid'
eliminate_rRNA = os.path.isfile(rnid_file)
del_rn={}
if eliminate_rRNA:
 print 'identified rRNA records will be removed, may add them back by editing rnid file'

for f in gff_r:
 added=[]
 candidates=sorted(gff_r[f],key = lambda x: x.qualifiers['score'],reverse=True)
 while len(candidates)>0:
  current = candidates.pop(0)
  if not any(overlapl(current.location,x.location) for x in added):
    added.append(current)
 rna_l = sorted(added,key = lambda x: x.qualifiers['score'])
 kingdom=''
 for each in rna_l:
  if each.type == 'rRNA':
    dom = each.qualifiers['product'].split('_')
    kingdom = dom[-1]
    each.qualifiers['product']=rna[dom[0]+dom[-1]]

    # check for incomplete rRNA here - if located at border, modify location.
    strandness = each.strand
    if each.location.start == 0:
     each.location = FeatureLocation(BeforePosition(0),each.location.end,strandness)
    if each.location.end == len_dict[f]: 
     each.location = FeatureLocation(each.location.start,AfterPosition(len_dict[f]),strandness)

  if each.type == 'CDS': #beware of this dirty hack if non-mt hmms are used in wise2!
    kingdom = 'mitochondria'
 rdict[f].features.extend(rna_l)
  
#it would be desirable to remove at least some overlapping features from records with DNA-level annotations.
# and edit/replace description accordingly
# cluster by overlap and retain the best scoring from each cluster? 
# Not so important anyway as these records will mostly be deleted anyway (as mitochondrial)
# but this is possible only after all features are added.

# consider simplifying re-annotation/filtering - actually it is not that useful
# what if all original annotations are retained? These would not be used in tbl2asn anyway...
# at least could use hierarchical clustering instead of (obsolete?) rclust...

 clustered_features = rclust(rdict[f].features) #this empties _all_ features, they have to be added back.
 for f_list in clustered_features:
  strong = 5
  for each in f_list:
   if each.type=='rRNA':
    strong=1
    each.qualifiers['GO_terms']=['0005840','0003735','0006412']
    #set of ribosomal rRNA GO terms (RFAM approved)
    # consider expanding the dbsearch to the whole RFAM and adding rfam2go dictionary support here instead 
    try:
      reference = each.qualifiers['db_xref']
      each.qualifiers['db_xref']=[reference]
    except:
      continue
  strong=min(len(f_list),strong) #add this number at most
  best=sorted(f_list, key = lambda x: x.qualifiers['score'],reverse=True)[0:strong]
  rdict[f].features.extend(best)

 cds_count=0
 longest_cds=0
 for anno in rdict[f].features:
  if anno.type=='CDS':
   cds_count +=1
   if len(anno) > longest_cds:
    longest_cds=len(anno)
 rdict[f].dbxrefs=[str(cds_count),str(longest_cds)] 

 if len(kingdom)>0:
   tid=k2i[kingdom]
   if kingdom in ['bacteria','archea']:
    rdict[f].description += ', bacterial contamination'
    del_dict[f] = True #is this necessary? Should be filtered out later anyway...
    #    a_tid=nx.ancestors(TAX,tid) 
    # these (two) sets could actually be pre-calculated just after reading-in the TAX database...
    # some 100x faster?
    # relatively small k2anc dictionary holding the sets is harmless
    a_tid=k2anc[kingdom]
    ctid=rdict[f].annotations.get('TaxID',None)
    
    # DEBUG: some contaminants slip through this filter
    # ctid is TaxID assigned by "proteins"
    # tid is TaxID from rDNA - if we are here its either 2 or archeal
    # a_tid is a set of taxons for a given kingdom (bacteria or archea) present in TAX (so quite large)
    # the procedure focuses on the situation where these two don't agree.
    # and ignores the case when they do agree. There is no followup.
    # these cases should be delt with here if needed. What if there was no TaxID assigned?
    # (ctid==None?). This should still work... None is not in any set...
    # For correct filtering, the data stored with the record should also be updated, not just itaxon
    # FILTERING IS BASED ON DATA STORED WITH THE RECORD!!!
    # so rdict[f].annotations['TaxID'] should be set to tid AND ['organism'] should be set to i2n[tid]
    
    if not (ctid in a_tid): 
     rdict[f].annotations['TaxID']=tid
     rdict[f].annotations['organism']=i2n[tid]
     try:
      itaxon[tid].add(f)
     except KeyError:
      itaxon[tid]=set([f])
     if ctid:
      itaxon[ctid].discard(f)
      if len(itaxon[ctid])==0:
       del itaxon[ctid]
   elif kingdom == 'mitochondria':
    rdict[f].description += ', mitochondrial transcript'
    if eliminate_mt:
     del_mt[f] = True
   else:
    a_tid=k2anc[kingdom]
    ctid=rdict[f].annotations.get('TaxID',None)
    # this is hardly conclusive... 
    # replacement of previously assigned TaxiD based on this is VERY problematic
    # consider NOT doing this at all.
    """
    if not (ctid in a_tid): 
     try:
      itaxon[tid].add(f)
     except KeyError:
      itaxon[tid]=set([f])
     if ctid:
      print 'eliminating', f, 'from taxid',i2n[ctid],'because rRNA points at',i2n[tid]
      itaxon[ctid].discard(f)
      if len(itaxon[ctid])==0:
       del itaxon[ctid]
    """
    #Not really conclusive... May cause inconsistency...
    # Could easily check if there is already taxon on the record and 
    # only attempt assigning this general one in certain case: in prokaryotic case does not matter, 
    # will always try (with good effect on filtering). 
    # But for Eukatyotic this may actually mask contaminants!!! Even if nothing is assigned (isoforms!).
    # There is room here for additional search...
    # RFAM with TAX data? SILVA? Or some nr subset?
    rdict[f].description += ', contains rRNA from '+ i2n[tid]
    # need additional assessment of true rDNA transcripts.
    # The only reliable way to do this is through external addition of rRNA records(s)
    # in a similar way as mitotranscripts are added.
    # the simple criterion is misleading and could only be used in the most obvious cases.
    # the simple criterion: three rRNA features in one record 
    # change: organism, description. Make sure it is not deleted, unless duplicated.
    # what about itaxon?? Should not change annotation and leave itaxon!
    # since no reliable source of taxonomy is preset in singular cm profiles, the only fair approach
    # would be to ignore this
    if eliminate_rRNA:
     del_rn[f]=True

    elif len(rdict[f].features) == 3 and cds_count == 0:
     rdict[f].description = 'rRNA transcript'
     rdict[f].annotations['organism'] = '.' # this is for counter-filtering
     deleted=del_dict.get(f,False)
     duplicated=dup_dict.get(f,False)
     if deleted and not duplicated:
      del del_dict[f]

#----real end of annotation here----#

t8=time.time()
#print t8-t6,'s'
print len(gff_r),'records modified by DNA data, filtering'

print ' reading abundance data'
# for multi-sample data sets there will be several different abundance estimates.
# data structure with dictionary of lists works and naturally extends to single sample case (one element list)
# The dictionary is keyed by transcript id and by sample name and stores TPM value( a float)
# add condition to quit if no abundance provided (too much depends on it now)

abu_t={}
abut_file=name[0]+'.fsa.abu.salmon.isoform.TMM.EXPR.matrix'
if not os.path.isfile(abut_file):
 abut_file=name[0]+'.abu.salmon.isoform.TMM.EXPR.matrix'
 if not os.path.isfile(abut_file):
   abut_file=name[0]+'-abundance-fsa-s/quant.sf'
   if not os.path.isfile(abut_file):
    abut_file=name[0]+'-abundance-s/quant.sf'
    if not os.path.isfile(abut_file):
     abut_file=None
if abut_file:
 print ' from:', abut_file
 with open(abut_file, 'r') as abu_h:
  samples = abu_h.readline()
  if 'matrix' in abut_file:
   sample_names = samples.strip().split('\t')
## some pre-processing would be desirable so that these names are shorter
   prefix_len=len(name[0])+1
   better_names=[]
   for each in sample_names:
    last_one=each.find('-abundance')
    better_names.append(each[prefix_len:last_one])
   sample_names=better_names
  else:
   sample_names = []  
  sample_names.append(name[0])
  for line in abu_h:
   abu = line.split('\t')
   key=abu[0]
   if 'matrix' in abut_file:
    value=[float(x) for x in abu[1:]]
    value.append(sum(value)/float(len(value)))
   else:
    value=[float(abu[3])]
   # HERE define the dictionary effectively changing list do dict
   # keys are in sample_names, values are in value
   # there should be a simple way to create dictionary from two lists...
   # yes, there is, with zip
   abu_t[key] = dict(zip(sample_names,value))
else:
 print "can't continue without abundance data"
 quit()
# Abundance tables should be calculated using fsa file for reliable final result
# However, run this abundance analysis after the final version of fsa is ready
# will use fsa-matrix, matrix, individual-fsa or individual (order of preference)

abu_d={}
abu_file=name[0]+'.fsa.abu.salmon.gene.TMM.EXPR.matrix'
if not os.path.isfile(abu_file):
 abu_file=name[0]+'.abu.salmon.gene.TMM.EXPR.matrix'
 if not os.path.isfile(abu_file):
  abu_file=name[0]+'-abundance-fsa-s/quant.sf.genes'
  if not os.path.isfile(abu_file):
   abu_file=name[0]+'-abundance-s/quant.sf.genes'
   if not os.path.isfile(abu_file):
    abu_file = None

if 'matrix' in abu_file and 'matrix' in abut_file:
 print ' samples:', len(sample_names)-1
elif 'matrix' in abu_file or 'matrix' in abut_file:
 print "inconsistent abundance files, can't continue."
 quit()

if abu_file:
 with open(abu_file, 'r') as abu_h:
  abu_h.readline()
  for line in abu_h:
   abu = line.split('\t')
   key=abu[0]
   if 'matrix' in abu_file:
    value=[float(x) for x in abu[1:]]
    value.append(sum(value)/float(len(value)))
   else:
    value=[float(abu[3])]
   abu_d[key] = dict(zip(sample_names,value))

print ' number of expressed genes:', len(abu_d)
print ' number of records with expression:',len(abu_t)
abu_taxon={}
for each in itaxon:
   if each in TAX: 
    gid_list=[]
    for rid in itaxon.get(each,set()):
     gid = '_'.join(rid.split('_')[:-1])
     gid_list.append(gid)
    abu=0
    for g_id in gid_list:
     abu+=abu_d.get(g_id,{name[0]:0})[name[0]] 
    if abu > 0: 
     abu_taxon[each]=abu
   else:
    print each,'not in taxonomy, check the taxonomy annotation routine or delete bad_taxon file'

taxon_set=set(list(itaxon)) 


if True:
# Loop over filtering if fsa was used, from here to tbl but was moved to later stage. 
#TO DO: remove bogus if and adjust indentation.


 print 'number of assigned taxons:',len(itaxon)

 bad_taxon=set()
 bad_taxon_file=name[0]+'.bad_taxon'
 if os.path.isfile(bad_taxon_file): 
 # once we are happy with bad_taxon file, there is no point in repeating the analysis of assigned taxons
  print ' reading list of potential contaminants'
  with open(bad_taxon_file, 'r') as fh:
   for line in fh:
    tid=line.strip().split('\t')[0].strip()
    if tid in TAX:
     bad_taxon.add(tid)
    else:
     if tid=='mgcode':
      mgcode = line.strip().split('\t')[1]
     else:
      print ' taxon unknown',line.strip()
     
 else:

  print 'analysing taxonomy of the sample'
  tt=time.time()
  ## limit TAX database here
  ## empty leaves and intermediate nodes will be eliminated
  node_list = []
  for x in TAX.nodes():
    if TAX.in_degree(x) == 0:
      node_list.append(x)
  node_list = set(node_list) - taxon_set 
  TAX.remove_nodes_from(node_list)
  while len(node_list) > 0:
   node_list = []
   for x in TAX.nodes():
    if TAX.in_degree(x) == 0:
     node_list.append(x)
   node_list = set(node_list) - taxon_set 
   TAX.remove_nodes_from(node_list)
#  print 'nodes in taxonomy after removal of empty leaves:', len(TAX)
  canon_list=['2','2157','10239','2698737','33090','2763','4751','9443','6231','50557','6656','7742','6157','35581','33208','2759']
  canonical_taxa=set(canon_list)
  node_list = []
  for x in TAX.nodes():
   if TAX.in_degree(x)==1 and TAX.out_degree(x)==1:
    node_list.append(x)
  node_list = set(node_list) - taxon_set - canonical_taxa
  while len(node_list) >0:
   remove_me=node_list.pop()
   pnodes=list(TAX.predecessors(remove_me))
   snodes=list(TAX.successors(remove_me))
   TAX.remove_nodes_from([remove_me])
   TAX.add_edges_from([(pnodes.pop(),snodes.pop())])
   node_list = []
   for x in TAX.nodes():
    if TAX.in_degree(x)==1 and TAX.out_degree(x)==1:
     node_list.append(x)
   node_list = set(node_list) - taxon_set - canonical_taxa
#  print 'nodes in taxonomy after removal of internal empty nodes:', len(TAX)

  ctaxon={}
  for tid in itaxon:
   clad= nx.descendants(TAX,tid) | nx.ancestors(TAX,tid) | set([tid])
   rid_list=itaxon.get(tid,set())
   for member in clad:
    curid=ctaxon.get(member,set())
    curid.update(rid_list)
    ctaxon[member]=curid

  print len(ctaxon),'taxons considered'
  tabu={}
  nabu={}
  for tid in ctaxon: # ? ctaxon or itaxon? After TAX reduction should not matter speed-wise but could influence whitelisting
   abu=[abu_t.get(x,{name[0]:0})[name[0]] for x in ctaxon[tid]]
   tabu[tid]=sum(abu)
   nabu[tid]=len(abu)
  scls=sorted(nabu, key=tabu.get, reverse=True)
  
  # sorting and decision should be based on the same stat!!
  # reporting should be limited to the top of the list... Consider hardcoding a limit scls length
  # also, saving this is pointless, onscreen reporting should be enough/better
  root_tid=scls[0]
  ttt=time.time()
  print ' taxons sorted by abu, saving clade-stats'#,ttt-tt,'s spent in taxonomy so far'
  # consider simplifying this csv
  #or eliminating, with fixed uniref taxids works flawlessly so further debugging is not needed
  nab=nabu[root_tid]
  tab=tabu[root_tid]
  fh=open(name[0]+'.clade-stats.csv','w')
  header=['taxon','tTPM', 'ttn','n%','a%','d']
  fh.write('\t'.join(header)+'\n') 
  for tid in scls:
   tbs=[i2n[tid],str(tabu[tid]),str(nabu[tid]),str(nabu[tid]/float(nab)),str(tabu[tid]/tab),str(nabu[tid]-nab)]
   fh.write('\t'.join(tbs)+'\n') 
  fh.close()

  # maxTPM=tabu[root_tid]*0.98
  #setting this too high will lead to low stringency, too low will remove legit relatives (too high stringency)
  #maxTno=nabu[root_tid]*0.955
  maxTno=tabu[root_tid]*0.98

  for tid in scls:
  # cuTPM=tabu[tid]
   cuTno=tabu[tid]
   if cuTno < maxTno:
    break
  #  else:
  ttid=tid

  # define bad_taxon database by subtracting whitelist from all taxons encountered 

  good_taxon = nx.descendants(TAX,ttid) | nx.ancestors(TAX,ttid) | set([ttid])
  bad_taxon = taxon_set - good_taxon
  print 'whitelisting the clade of:',i2n[ttid],cuTno/float(tabu[root_tid]),'(implied purity), check clade file for details'
  # propose a bad_taxon file augmented with abundance data (already precalculated).
  del tabu, nabu, ctaxon

  ## It would be nice to guess the mt genetic code based on this. It is in nodes.dmp, in sixth (?) column, with tid at column 0.
  ## scan through the file and store the value
  with open('nodes.dmp','r') as nh:
   for ln in nh:
    col=ln.split('\t|\t')
    if col[0]==ttid:
     mgcode=col[8]
     break
  print 'will use',mgcode,'table for mitochondrial genetic code'
  # canonical taxa are defined for augmenting bad_taxon. 
  # May add more subgroups to Metazoa (33208) but it works pretty well 
  # consider guessing true contaminants and false positives automatically,
  # by splitting the list of canonical taxa into "suspicious" and "benign"
  # (low priority, manual sorting out of the augmented bad_taxon file is quite fast)

  tax_h = open(bad_taxon_file,'w')
  taxon_l=sorted(bad_taxon, key=abu_taxon.get, reverse=True)
  for each in taxon_l:
    dsc=nx.descendants(TAX,each)
    aug_set = dsc & canonical_taxa
    tax = i2n[each]
    for tid in canon_list:
     if tid in aug_set:
      parent=i2n[tid]
      tax = ' '.join([tax,parent])
      break
    abu = abu_taxon.get(each,0)
    tax_h.write('\t'.join([each,tax,str(abu)])+'\n')
  tax_h.write('mgcode\t'+mgcode+'\n')
  tax_h.close()
 
# end of classification here - back to taxon db filtering
# - bad_taxon is a set of taxon names, either generated or read-in
 sum_bad_abu=0
 #Now can filter (=add to del_dict) 
 # records from itaxon which are in the bad_taxon iterable. Isoform check will follow.
 for each in bad_taxon:
   abu = abu_taxon.get(each,0)
   sum_bad_abu+=abu
   r_list=itaxon.get(each,[]) 
   # get() needed for experimentation on already filtered data sets, harmless otherwise
   for r_id in r_list:
     del_dict[r_id]=True
     rdict[r_id].description += ' suspicious taxon'

 print len(bad_taxon), 'bad taxons cover approx.', sum_bad_abu//10000,'% of expression, examine and adjust bad_taxon file'
 del abu_taxon, itaxon, TAX, i2n

 print 'filtering'

 print ' gathering deletion candidates:', len(del_dict)

 # it is better not to delete all isoforms blindly. Mitochondrial are OK to delete. 
 # Enforced or duplicate are NOT to be deleted and should not be in del_dict yet.
 # suspicious bad_taxons may be deleted ONLY if none of isoforms is positively good and some are bad.
 # this must be checked. 
 #Final assessment of the overall contamination (TPM) would be nice but is not mandatory
 # (not implemented yet)
 
 # iterate over rdict and build dictionary of isoforms, keyed by gid, with list of rids (gdict)
 # identify problematic lists of rids - having the mentioned condition.
 if eliminate_mt:
  with open(id_file,'r') as id_h:
   line= id_h.readline()
   id_n=line.split('\t')
   mitogid='_'.join(id_n[1].strip().split('_')[:-1])
 else:
  mitogid=''
 if eliminate_rRNA:
  with open(rnid_file,'r') as id_h:
   line= id_h.readline()
   id_n=line.split('\t')
   rngid='_'.join(id_n[1].strip().split('_')[:-1])
 else:
  rngid=''
 del_short={}
 gdict={}
 #final filtering is based on rdict contents, not itaxon nor del_dict !!! Beware !!!
 # rdict must be consistent with itaxon!
 for f in rdict:
  if len(rdict[f]) < 200:
    rdict[f].description += ' too short'
    del_short[f] = True
  g_id = '_'.join(f.split('_')[:-1])
  try:
   gdict[g_id].append(f)
  except KeyError:
   gdict[g_id]=[f]
 for gid in gdict:
  #build a list of "classifiers": good/bad and taxon/not taxon for each gene (gid)
  disqualifiers=[del_dict.get(x,False) for x in gdict[gid]]
  mito=[del_mt.get(x,False) for x in gdict[gid]]
  rrna=[del_rn.get(x,False) for x in gdict[gid]]
  if any(disqualifiers):
   # classifiers should better be derived from itaxon...? Orelse strive to keep these consistent!!!
   classifiers=[rdict[x].annotations.get('TaxID',None) for x in gdict[gid]]
   good_classifiers=[x for x in classifiers if x]
   firm_disqualifiers=[x for x in disqualifiers if x]
   if len(good_classifiers) == len(firm_disqualifiers):
    for rid in gdict[gid]:
     if not del_dict.get(rid, False):
      del_dict[rid]=True
      rdict[rid].description += ' suspicious gene'
   elif any(mito) or gid == mitogid:
    for rid in gdict[gid]:
     if not del_dict.get(rid, False):
      del_dict[rid]=True
      rdict[rid].description += ' mitotranscript'
   elif any(rrna) or gid == rngid:
    for rid in gdict[gid]:
     if not del_dict.get(rid, False):
      del_dict[rid]=True
      rdict[rid].description += ' rRNA'
   else:
    for rid in gdict[gid]:
     if del_dict.get(rid, False):
      del del_dict[rid]
 del_dict.update(del_mt)
 del_dict.update(del_short)
 del_dict.update(del_rn)
 print ' check isoforms and eliminate very short sequences:', len(del_dict)      

 del_dict.update(dup_dict) # after isoform, before removal 
 if len(dup_dict) > 0:
  print ' with duplicates:', len(del_dict)   

 # move this to the end to prevent unintentional deleting of isoforms
 del_file=name[0]+'.delete'
 if os.path.isfile(del_file):
  with open(del_file,'r') as del_h:
    for line in del_h:
     transcript = line.split('\t')[0].strip()
     del_dict[transcript]=True
     rdict[transcript].description += ' enforced deletion'
  print ' with deletions enforced by file:', len(del_dict)

 d_list=sorted(del_dict, key = lambda x: abu_t.get(x,{name[0]:0})[name[0]], reverse=True)
 eliminated=open(name[0]+'.eliminated', 'w')
 for f in d_list:
    SeqIO.write(rdict[f],eliminated,'gb')
 eliminated.close()
 if not ('fsa' in abu_file):
  # store detailed contamination report for all samples.
  #only possible for raw abundaance data
  empty_dict=dict(zip(sample_names,[0]*len(sample_names))) #for sorting and adding
  header='\t'.join(sample_names)
  with open(name[0]+'.eliminated.csv','w') as csv:
   csv.write('transcript\tproduct\ttaxon\tlen'+header+'\n')
   for f in d_list:
    r=rdict[f]
    abu=abu_t.get(f,empty_dict)  # returns dictionary
    abu='\t'.join([str(abu[x]) for x in sample_names])
    line_list=[r.id,r.description,r.annotations.get('organism',''),str(len(r)),abu]
    line = '\t'.join(line_list)+'\n'
    csv.write(line)

 if len(del_dict) > 0:
  print len(del_dict), 'eliminated records saved in gb file'
 del d_list
 # THIS IS INDENTED only because of bogus if, to be removed.

print 'Filtering by abundance'
## prefiltering based on matrix of gene expressions: remove genes and isoforms
## which are not expressed
# =>Directly remove all records not mentioned in abundance table as alternative to the following.
if not ('fsa' in abu_file):
 print ' based on raw expression matrix: all stats are provisional'
 remain_candidates=set(rdict)-set(del_dict)
 sample_set={}
 for sample in sample_names:
  working_list=[]
  for f in remain_candidates:
   g_id = '_'.join(f.split('_')[:-1])
   length=len(rdict[f])
   try:
    orfl=rdict[f].dbxrefs[1]
   except:
    orfl=0
 #  if abu_d[g_id][sample] > 1 and abu_t[f][sample] >0.1:
 #  if abu_t[f][sample] >10 or (abu_d[g_id][name[0]]<10 and abu_t[f][sample] > 1 and length > 500 and orfl > 200):
 #  conversion to get() was necessary to process re-filter after fsa abu.
 #  BUT this should never be needed. If filtering, then only based on complete set - there should be no KeyErrors!
 # reverse, this is impossible/impractical to keep. TO allow mismatch between abu and rdict, return to get()
 #  if abu_t[f][sample] >= 10 or (abu_d[g_id][name[0]]<20 and abu_t[f][sample] > 1 and length > 400 and orfl > 150):
 #  if abu_t[f][sample] > 10 or (abu_d[g_id][name[0]]<10 and abu_t[f][sample] > 1 and length > 500 and orfl > 200):
   
   if abu_t.get(f,{sample:0})[sample] >=10 or (abu_d.get(g_id,{name[0]:0})[name[0]]<20 and abu_t.get(f,{sample:0})[sample] > 1 and length > 400 and orfl > 150):
 #  if abu_t.get(f,{sample:0})[sample] >1:# or (abu_d.get(g_id,{name[0]:0})[name[0]]<10 and abu_t.get(f,{sample:0})[sample] > 1 and length > 500 and orfl > 200):

   ##### THIS IS THE KEY LINE OF CODE### TUNE the values up for more aggressive filtering
   ##### think and experiment for much smaller final db. Check results with BUSCO
   ##### transrate should jump significantly with each increase of stringency

   # current setting: is obviously expressed or has nice structure
   # there will be some obviously expressed with bad structure though... 
   # maybe for those obvious, if there are annotations, these should favor nice structures?
   # No way to be sure there is an isoform with nice structure, though...
   # Essentially, the difference between high-expression and low-expression should be that for high expression
   # records without annotations should pass, for low expression they should not pass.
   # But if there are annotations, they should be good enough.

    working_list.append(f)
  sample_set[sample]=set(working_list)
 # To exclude the average sample just remove the name[0] set...
 # and iterate over sample sets to get union and del sets

 if len(sample_set) >1:
  del sample_set[name[0]] 
  # sample set is now a dictionary with only real samples in it and so if average is to be used in filtering
  # it must be included in the condition above. In any case there must be some expression of anything kept
  # in at least one sample so the average "sample" is not necessary at this step.
 union_set=set()
 for sset in sample_set:
  union_set = union_set.union(sample_set[sset])
 del_set = remain_candidates - union_set

 del sample_set, remain_candidates, union_set

## denovo abundance filtering ends here; alternative filtering based on the contents of abu_t follows

else:
 print ' based on secondary expression matrix'
# del_dict={}
 kept_records=set(list(abu_t))
 all_records=set(list(rdict))
 del_set=all_records - kept_records

# end of gathering deletion candidates
# execute deletion

eliminated=open(name[0]+'.poor', 'w')

for f in del_set:
 del_dict[f]=True
 #save this for troubleshooting
 SeqIO.write(rdict[f],eliminated,'gb')
eliminated.close()
print len(del_set),' poorly expressed records saved in poor file'
 
for f in del_dict:
 del rdict[f]
del del_dict, del_set #may add some more but it is not critical.

#no more deletions, don't need all these sets anymore.
# replacement of mitochondrial records must still be completed

# small loop to optionally add back rRNA record(s)
# no need(?) for special treatment, manual control of gbk file contents is preferable
if eliminate_rRNA:
 with open(rnid_file,'r') as id_h:
  for line in id_h:
   id_n=line.split('\t')
   gbk_name = id_n[0] + '.gbk'
   print 'adding back rRNA data from',gbk_name
   new_rec_h = SeqIO.parse(gbk_name,"genbank")
   r = new_rec_h.next()
   r.id = id_n[1].strip()
   for f in r.features:
    for qual in f.qualifiers:
     if qual=='product':
      f.qualifiers[qual] = f.qualifiers[qual][0]
     elif qual == 'score':
      f.qualifiers[qual] = float(f.qualifiers[qual][0])
   rdict[r.id] = r

# everything below this point is for reporting of the final, filtered or pre-filtered data set, 
# with the assumption that mitotranscripts will fit abu more or less OK.

print len(rdict), 'records left after all filtering steps'
print 'preparing tbl & fsa for submission'

# accommodate dr name suggestions here. Optionally (file exists) preload edited dr suggestions 
# and use it as a dictionary to correct product names (only) in records/features.
# keys must be composed of record.id and feature location.
# values should be calculated with real data at hand and replaced if present in dictionary (d.get(value,unchanged))

sugfile=name[0]+'.sug'
sug_dict={}
if os.path.isfile(sugfile):
 print 'acquiring name suggestions:',
 with open(sugfile,'r') as sug_h:
    for line in sug_h:
     if len(line) <10:
      break
     sug=line.split('\t')
     sug_dict[sug[2]]='hypothetical protein'
 print len(sug_dict),'from preliminary tbl2asn'

# consider augmenting, not discarding entirely (?) orelse using alternative files to avoid redundant IO
edited_tbl_file=name[0]+'.tbl_edited'
if os.path.isfile(edited_tbl_file):
 sug_dict=parse_tbl(edited_tbl_file)
 print len(sug_dict),'suggestions overwritten by edited tbl' 

allowed_tags=['product','db_xref']
allowed_types=['CDS','rRNA']
optional_db=['PANTHER']

tbl=open(name[0]+'.tbl','w')
fsa=open(name[0]+'.fsa','w')

final_annotated=0

for record in rdict:
    r=rdict[record]
    # reversed CDS are not welcome in TSA. 
    # Therefore records with best annotations in minus frame must be reversed before fsa and tbl
    ### All remaining features in minus will not be annotated in tbl ###
    try:
     strandness=sorted(r.features,key = lambda x: float(x.qualifiers['score']), reverse=True)[0].strand
    except:
     strandness = 1
    if strandness == -1:
      r = r.reverse_complement(id=r.id,description=r.description)
    dscript = r.description
    r.description = '[moltype=transcribed_RNA] [tech=TSA]'
    SeqIO.write(r,fsa,'fasta')
    r.description = dscript
    if len(r.features) > 0:
      header = '>Feature lcl|' + str(r.id) +'\n'
      contents=''
      for f in r.features:
       if f.type in allowed_types and float(f.qualifiers['score']) > 20 and not (f.type == 'CDS' and len(f) <150):
        if f.strand == 1:
         contents += str(f.location.start+1)+'\t'+str(f.location.end)+'\t'+str(f.type)+'\n'
         loc_str= str(f.location.start+1)+'-'+str(f.location.end)
         for key,val in f.qualifiers.iteritems():
          if key in allowed_tags:
           if key=='product':
            better_val = sug_dict.get(r.id+':'+loc_str, val)
            val = better_val
            if val.upper() == val:
             val+=' protein'
            val = [val]
           for each in val:
            db_id=each.split(':')[0]
            if db_id in optional_db:
             contents+='\t\t\tnote\tdb_xref_other='+str(each)+'\n'
            else: 
             contents+='\t\t\t'+str(key)+'\t'+str(each)+'\n'
      if contents:
       tbl.write(header+contents) 
       final_annotated+=1

t10=time.time()
# print ' annotated:',final_annotated
print t10-t8,'s\n checking for mito'
mt_set=set()
if eliminate_mt: # this condition should be simplified (?) or more precise. Filtered data do not imply id-mapping... But this file could also be empty.
  non_list_single_terms=['structure','anticodon','product']
  id_h=open(id_file,'r')
  for line in id_h:
   id_n=line.split('\t')
   gbk_name = id_n[0] + '.gbk'
   print 'adding data from',gbk_name
   new_rec_h = SeqIO.parse(gbk_name,"genbank")
   r = new_rec_h.next()
   sequence_len = len(r)
   r.id = id_n[1].strip()
   rdict[r.id] = r 
   ## may need to secure this to avoid unintentional overwriting 
   ## the problem is that it is too late to update fsa - some records may get duplicated!
   ## Need to store these id in a mt_set for normalization
   mt_set.add(r.id)
   final_annotated+=1
   for f in r.features:
    for qual in f.qualifiers:
     if qual in non_list_single_terms:
      f.qualifiers[qual] = f.qualifiers[qual][0]
    if f.type == 'tRNA':
     try:
      f.qualifiers['structure'] = f.qualifiers['structure'].replace(' ','')
     except:
      pass
   # initial correction of anticodon location is based on structure so do need structure
   # reversal of anticodon is relying on internal Biopython reversal of features (for trn on minus strand)
   for f in r.features:
    if f.type == 'tRNA':
       parts = f.qualifiers['anticodon'][1:-1].split(',',1)
       offset = f.qualifiers['structure'].index('___')
       new_loc = FeatureLocation(f.location.start + offset,f.location.start+offset+3, f.strand)
       f.qualifiers['anticodon'] = '(pos:'+ loc2gb(new_loc)+','+parts[1]+')'
    elif f.type== 'CDS':
       f.qualifiers['transl_table']=mgcode
 
   allowed_tags=['product', 'anticodon','transl_table']
   allowed_types=['CDS','rRNA','tRNA','polyA_site']
   tbl.write('>Feature lcl|' + r.id+'\n')
   try:
     strandness=sorted(r.features,key = lambda x: len(x), reverse=True)[0].strand
   except:
     strandness = 1
   if strandness == -1:
    fake_anticodon=[]
    for f in r.features:
     if f.type=='tRNA':
      parts= f.qualifiers['anticodon'][1:-1].split(',',1)
      orig_loc = gb2loc(parts[0])
      fake_anticodon.append(SeqFeature.SeqFeature(orig_loc,'misc',{}))
    r.features.extend(fake_anticodon)
    r=r.reverse_complement(id=r.id,description=r.description)
    for f in r.features:
     if f.type=='tRNA':
      parts= f.qualifiers['anticodon'][1:-1].split(',',1)
      for g in r.features:
       if g.type=='misc' and g.location.start in f.location:
        new_loc=g.location
        break
      f.qualifiers['anticodon']='(pos:'+loc2gb(new_loc)+','+parts[1]+')'
      # record reversal with anticodon recalculation.
      # first decode location of anticodon
      # then define new feature in this location (save the rest of the anticodon string)
      # now revcomplement the record. Get the location of the previously defined feature
      # modify the anticodon with the new location and the rest of the string.
      # deleting of the temporary feature is not needed
   for f in r.features:
    if f.type in allowed_types:
     exception = None
     fx = f.location
     phase = len(fx) % 3
     if f.strand == -1:
      print 'features on reverse are ignored' # this should not happen, inspect gbk (defining mt pieces) files if reported
     else:
      if f.type == 'CDS' and phase <> 0:
        if phase == 1:
          last_codon = FeatureLocation(fx.parts[-1].end-1,fx.parts[-1].end, fx.strand)   
        elif phase == 2:
          fx = diminute_end(f.location)                                                        ## Workaround tbl2asn bug: shorten the CDS by one
          last_codon = FeatureLocation(fx.parts[-1].end-1,fx.parts[-1].end, fx.strand)         ## should be end-1, if not for the bug
        exception = '\t\t\ttransl_except\t(pos:' + loc2gb(last_codon) + ",aa:TERM)\n\t\t\tnote\tTAA stop codon is completed by addition of 3' A residues to mRNA\n"
      tbl.write(str(fx.parts[0].start+1)+'\t'+str(fx.parts[0].end)+'\t'+str(f.type)+'\n')
      if len(fx.parts) > 1:
       for i in range(1,len(fx.parts)):
           tbl.write(str(fx.parts[i].start+1)+'\t'+str(fx.parts[i].end)+'\n')
     for key,val in f.qualifiers.iteritems():
      if key in allowed_tags:
       tbl.write('\t\t\t'+str(key)+'\t'+str(val)+'\n')
     if exception:
      tbl.write(exception)
   dscrpt=r.description
   r.description = '[moltype=transcribed_RNA] [tech=TSA] [location=mitochondrion] [mgcode='+mgcode+']'
   SeqIO.write(r, fsa, 'fasta')
   r.description = dscrpt
  id_h.close()

print final_annotated,'records annotated, closing fsa & tbl files'
print ' preparing csv and GenBank with summary'
tbl.close
fsa.close

# tbl and fsa data are saved.
# the mito records are also in via id-file mapping
# remember, to have mitotranscripts under one "gene" they must differ in "isoform" only, not "gene"
# so use 'TRINITY_DNxxx_c0_g0_in style, not 'TRINITY_Dxxx_c0_gn_i1' in records listed in id_file.

# prepare csv file with record_id, description and organism
# include additional data not available in clc: 
# number of annotations, number of annotated  CDS (orfs) as well as length
# augmented with abundance information from the quant/matrix file and GO terms
# in this loop populate also term_dict to be used in GO reporting later (index encountered terms in rdict, to iterate over terms)
# Actually, it would be desirable to move term_dict to the first loop. Or merge the loops again (they were on once)
# The reason to split was the issue with samples which is no longer there (IIRC)

csv=open(name[0]+'.csv','w')
TPM='\t'.join(sample_names)
line = '\t'.join(['transcript ID','likely product',TPM,'taxon','length','longest','#ORF','#annot','GO'])+'\n'
csv.write(line)

empty_dict=dict(zip(sample_names,[0]*len(sample_names))) #for sorting and adding
very_empty_dict=dict(zip(sample_names,[None]*len(sample_names))) #for counting
r_list=sorted(rdict, key = lambda x: abu_t.get(x,{name[0]:0})[name[0]], reverse=True)
r_list=sorted(r_list, key = lambda x: abu_d.get('_'.join(x.split('_')[:-1]), {name[0]:0} )[name[0]], reverse=True)
term_dict={} # dictionary keyed by represented terms, to keep sets of records, after expansion 
for record in r_list:
    r=rdict[record]
    abu=abu_t.get(record,empty_dict)  # returns dictionary
    abu='\t'.join([str(abu[x]) for x in sample_names]) # x is a key better get it from sample_list than from abu (order!)
    try:
     cds_count=r.dbxrefs[0]
     longest_orf=r.dbxrefs[1]
    except:
     cds_count='0'
     longest_orf='0'
    dscript = r.description
    GO_set=set()
    for f in r.features:
     if f.type=='CDS':
      try:
       GO_set.update(f.qualifiers['GO_terms'])
      except:
       pass
    tset=set()
    for each in GO_set:
     dsc=nx.descendants(rGO,each)
     tset.update(dsc)
    GO_set.update(tset) 
    # this one is used in GO reporting (not to repeat nx expansion) via term_dict.
    for term in GO_set:
     try:
      term_dict[term].add(record)
     except KeyError:
      term_dict[term]=set([record])
    GO_terms =','.join(list(GO_set - set(BAD)))# This one should be filtered for non-annotation database
    line_list=[r.id,dscript,abu,r.annotations.get('organism',''),str(len(r)),longest_orf,cds_count,str(len(r.features)),GO_terms]
    line = '\t'.join(line_list)+'\n'
    csv.write(line)
csv.close
print 'record list sorted by raw average gene-level abundance saved in csv'
print 'saving complete records in GenBank file'
t=time.time()
output_name=name[0]+'.gb'
fh = open(output_name,'w')
for x in r_list: # as a bonus output these in the same order, sorted by average abundance!
  record=rdict[x] 
  SeqIO.write(record,fh,'gb')
fh.close()
tt=time.time()
print tt-t,'s'
print len(r_list), 'records saved:', output_name
t10=time.time()
print t10-t1,'s so far'

## Dominance of GO landscape by  mitotranscripts is unjustified but real.
#  Data saved in csv are affected but can be manipulated downstream; can't update those easily since term_dict is build in the same loop as csv.
#  For Venn & GO reporting there is a possibility to adjust mt transcript expression so that it is better justified.
#  Two possibilities:
#  -exclude mitochondrial records from term_dict. Relatively easy - just remember mt_set and do the filtering at each term, using set math
#  OR:
#  -normalize abu_t and abu_d for mitotranscripts and mitogene (one gid is currently used by convention)
#   +collect set of oxphos gene terms 0005747, 0005749, 0005750, 0005751, 0016470 or single inner mt membrane term 0098800 from term_dict
#   +calculate nuclear subset by subtracting mt set
#   +calculate mitochondrial subset by intersecting with mt set
#   +get list of abu_t for mt and n sets
#   +calculate average TPM for nuclear and for mt set
#   +divide nuclear by mt to get scaling factor
#   +multiply abu_t by scaling factor for each in mt set (not only intersection, these will affect rRNA primarily, an intentional side effect) 
#   +multiply abu_d of mt gid by the same scaling factor
if len(mt_set) > 0 and len(abu_t) == len(rdict):
 print ' normalizing mitochondrial expression for GO terms'

# inner=term_dict.get('0098800',set())
 complexes=['0005747', '0005749', '0005750', '0005751', '0016470']
 OXPHOS=set()
 for c in complexes:
  OXPHOS = OXPHOS.union(term_dict.get(c,set()))
# n_inner = inner - mt_set
 n_OXPHOS = OXPHOS - mt_set
# mt_inner = inner & mt_set
 mt_OXPHOS = OXPHOS & mt_set
# print 'number of inner mt membrane transcripts:', len(inner),len(n_inner),len(mt_inner)
# print 'number of OXPHOS transcripts:', len(OXPHOS),len(n_OXPHOS),len(mt_OXPHOS)
# mt_OX_abu=[abu_t.get(x,empty_dict)[name[0]] for x in mt_OXPHOS if abu_t.get(x,empty_dict)[name[0]]>0]
# mt_inner_abu=[abu_t.get(x,empty_dict)[name[0]] for x in mt_inner if abu_t.get(x,empty_dict)[name[0]]>0]
# n_OX_abu=[abu_t.get(x,empty_dict)[name[0]] for x in n_OXPHOS if abu_t.get(x,empty_dict)[name[0]]>0]
# n_inner_abu=[abu_t.get(x,empty_dict)[name[0]] for x in n_inner if abu_t.get(x,empty_dict)[name[0]]>0]
# print 'TPM covered by and number of and median:'
# print 'inner mt membrane transcripts (nuclear,mt):', sum(n_inner_abu),sum(mt_inner_abu),len(n_inner_abu),len(mt_inner_abu),statistics.median(n_inner_abu),statistics.median(mt_inner_abu)
# print 'OXPHOS transcripts (nuclear,mt):', sum(n_OX_abu),sum(mt_OX_abu),len(n_OX_abu),len(mt_OX_abu),statistics.median(n_OX_abu),statistics.median(mt_OX_abu)
# print 'scaling (inner, OX):', (sum(mt_inner_abu)*len(n_inner_abu))/(sum(n_inner_abu)*len(mt_inner_abu)), (sum(mt_OX_abu)*len(n_OX_abu))/(sum(n_OX_abu)*len(mt_OX_abu))
# print 'median-based scaling', statistics.median(mt_inner_abu)/statistics.median(n_inner_abu),statistics.median(mt_OX_abu)/statistics.median(n_OX_abu)
 # check also at gene level... So convert sets of tids into sets of gids...
# n_inner_g=set()
# for tid in n_inner:
 # gid='_'.join(tid.split('_')[:-1])
#  n_inner_g.add(gid)
# mt_inner_g=set()
# for tid in mt_inner:
#  gid='_'.join(tid.split('_')[:-1])
#  mt_inner_g.add(gid)
 n_OXPHOS_g=set()
 for tid in n_OXPHOS:
  gid='_'.join(tid.split('_')[:-1])
  n_OXPHOS_g.add(gid)
 mt_OXPHOS_g=set()
 for tid in mt_OXPHOS:
  gid='_'.join(tid.split('_')[:-1])
  mt_OXPHOS_g.add(gid)
 mt_OX_abug=[abu_d.get(x,empty_dict)[name[0]] for x in mt_OXPHOS_g if abu_d.get(x,empty_dict)[name[0]]>0]
# mt_inner_abug=[abu_d.get(x,empty_dict)[name[0]] for x in mt_inner_g if abu_d.get(x,empty_dict)[name[0]]>0]
 n_OX_abug=[abu_d.get(x,empty_dict)[name[0]] for x in n_OXPHOS_g if abu_d.get(x,empty_dict)[name[0]]>0]
# n_inner_abug=[abu_d.get(x,empty_dict)[name[0]] for x in n_inner_g if abu_d.get(x,empty_dict)[name[0]]>0]
# print 'TPM covered by and number of and median:'
# print 'inner mt membrane genes (nuclear,mt):', sum(n_inner_abug),sum(mt_inner_abug),len(n_inner_abug),len(mt_inner_abug),statistics.median(n_inner_abug), statistics.median(mt_inner_abug)
# print 'OXPHOS genes (nuclear,mt):', sum(n_OX_abug),sum(mt_OX_abug),len(n_OX_abug),len(mt_OX_abug),statistics.median(mt_OX_abug),statistics.median(n_OX_abug)
# print 'gscaling (inner, OX):', (sum(mt_inner_abug)*len(n_inner_abug))/(sum(n_inner_abug)*len(mt_inner_abug)), (sum(mt_OX_abug)*len(n_OX_abug))/(sum(n_OX_abug)*len(mt_OX_abug))
# print 'median-based scaling', statistics.median(mt_inner_abug)/statistics.median(n_inner_abug),statistics.median(mt_OX_abug)/statistics.median(n_OX_abug)
 # gene-based median of OXPHOS is the best.
 #scaling_factor=statistics.median(mt_OX_abug)/statistics.median(n_OX_abug)
 #print ' average mitotranscript scaling factor:', scaling_factor
 """
 for s in sample_names:
  for tid in mt_set:
   abu_t[tid][s]=abu_t[tid][s]/scaling_factor
  gid='_'.join(tid.split('_')[:-1])
  abu_d[gid][s]=abu_d[gid][s]/scaling_factor
 """  
print len(term_dict),'GO terms with expression'

# iterate over samples but always use the whole rdict 
# the differentiating part is the abu dictionary second key (sname)
# iterate through sample_names (to keep the same order)
lvl={}
for sname in sample_names:
 if len(abu_t) == len(rdict): # better compare the sets but check if this works at all (and if precision is worth the effort)
  mt_OX_abug=[abu_d.get(x,empty_dict)[sname] for x in mt_OXPHOS_g if abu_d.get(x,empty_dict)[sname]>0]
  n_OX_abug=[abu_d.get(x,empty_dict)[sname] for x in n_OXPHOS_g if abu_d.get(x,empty_dict)[sname]>0]
  scaling_factor=statistics.median(mt_OX_abug)/statistics.median(n_OX_abug)
  print ' mitotranscript scaling factor:', scaling_factor
  #print 'ignoring mitoscaling'
  
  for tid in mt_set:
   abu_t[tid][sname]=abu_t[tid][sname]/scaling_factor
  gid='_'.join(tid.split('_')[:-1])
  abu_d[gid][sname]=abu_d[gid][sname]/scaling_factor
  
 else:
   print 'provisional ',
 print 'GO for',sname
 # term_dict and per-sample abu_t and abu_d lists for each term in set
 # are used to calculate  sums of numbers for all encountered terms
 # all terms are reported except the (presumabley three) root terms
 # each under the respective 0-level term name

 for j in range(len(class_name)): 
 # j iterates over (presumabely four: root 0 and each of root categories 1 2 3) GO sets of terms
 # uses term_dict keyed by expanded term holding sets of directly and indirectly assigned tids
 # filter these sets by sample-specific abundance (use abu_t, not abu_d in initial filtering)
 # to get sample-specific GO sets and stats
  cname=class_name[j]
# for cname in classify_d['root']:
  all_terms=classify_d[cname]
  if sname==name[0]:
    root_fname=name[0]
  else:
    root_fname=name[0]+'.'+sname

  sorted_terms = sorted(all_terms, key = lambda x: len(term_dict.get(x,[])), reverse=True)
  # sorting allows limiting iteration to terms encountered in the whole data set.  
  # but this sorting is not sample-specific!
  # consider storing data for each sample in memory and sorting them individually before saving.
  # would also allow cleaner handling of opened files in the context of j...
  # actually j==0 set (Venn) is only needed for the average sample, for each tissue its kind of useless
  # and likewise j>0 sets are not needed for the average sample unless its the only sample... 
  # so it would be better to first iterate over j (or directly class_name, with 'root' cname equivalent to j==0), then over samples.
  # TO DO figure how to ensure the 'root' class_name is processed last, without using numerical ids
  gdensity = {}
#  lvl={} # re-calculating levels for each sample seems excessive. However this is done only for represented terms, in theory each sample could have different ones
         # consider hybrid solution - with try or with get() checking for existencce before nx and without resetting.
  gabu={}
  tabu={}
  glst={}
  # get the records from term_dict.get()
  # and then filter the list by abu_t.get(record)[sample])  
  for i in range(len(sorted_terms)): # need i to save data for Venn (j=0 & i=0,1,2). Otherwise term IDs would have to be used (also doable).
   term=sorted_terms[i]
   id_set=term_dict.get(term,set())
   tid_set=set(x for x in id_set if abu_t.get(x,empty_dict)[sname] > 0)
   transcript_number = len(tid_set)
   if transcript_number == 0:
    break # end iteration at first empty term (every subsequent term will also be empty)
   else:
    #the tricky part: convert list of records to list of abundances, 
    #the list should not include records not expressed in a given sample
    #moreover, do it (also?) at gene level, so first construct set of g_ids (!redundancy)
    # going for gene level is problematic, don't rely on it, it may blurr the differences between tissues.
    # explore both options and compare the results. 
    ### Preferably store also the complete term_dict in some form.
    if j >0:     # don't store root stats in csv, generate Venn diagrams instead, for root just store tid sets.
     gid_set=set()
     for each in tid_set:
      g_id = '_'.join(each.split('_')[:-1])
      gid_set.add(g_id)
     # dont iterate, build the list directly
     gabu[term]=[abu_d[x][sname] for x in gid_set]
     tabu[term]=[abu_t[x][sname] for x in tid_set]
     glst[term]=sorted(gid_set, key = lambda x: abu_d.get(x,empty_dict)[sname], reverse=True)[:min(len(gid_set),4)]
     # gene is guaranteed to have expression if transcript has it so no additional conditions
     #calculate root distance. we do have current term, its root is in the class_name.
     # so the distance in rGO between class_name and term is needed. Start at term, end at class_name.
     # to avoid redundancy check if this was already calulated and calculate only if not
     if not lvl.get(term,False):
      lvl[term] = nx.shortest_path_length(rGO,term,cname)
     gdensity[term]=sum(gabu[term])//len(gabu[term])
#     gdensity[term]=statistics.median(gabu[term])
    else: #store sets of ids for root (for Venn diagrams) but really need to do this for last sample only
     if i == 0:
      A=tid_set
     elif i == 1:
      B=tid_set
     elif i == 2:
      C=tid_set
  if j > 0 and (sname!=name[0] or len(sample_names)==1):
#  if sname!=name[0] or len(sample_names)==1:
   #don't report csv for the last sample, unless it is the only sample
   #alternatively report csv also for root (j==0), but also skip the average sample
   fh = open(root_fname+'.GO_'+ref[cname].strip()+'.csv','w')
   fh.write('term\tlvl\ttr num\tsum tabu\tgene num\tsum gabu\tgdensity\ttranscripts\n')
   sample_terms=sorted(gdensity, key=gdensity.get, reverse=True)
   for term in sample_terms:
    tskpts=','.join(glst[term])
    tbw=[ref[term],str(lvl[term]),str(len(tabu[term])),str(sum(tabu[term])),str(len(gabu[term])),str(sum(gabu[term])),str(gdensity[term]),tskpts]
    fh.write('\t'.join(tbw)+'\n')
   fh.close()
   #pick best sample_terms for display. Three best(?) of levels 1-4. 
   #Easy to extend to more levels or more ranks. This will go for all root terms (so x3)
   dt1=[x for x in sample_terms if lvl[x]==1][0]
   dt2=[x for x in sample_terms if lvl[x]==2][0]
   dt3=[x for x in sample_terms if lvl[x]==3][0]
   dt4=[x for x in sample_terms if lvl[x]==4][0]
   if j in [1,3]:# or sort and display best/threshold?
    densest_terms=[dt2,dt3,dt4]
   else:
    densest_terms=[dt1,dt2,dt3]
   for term in densest_terms:
    print '\t',ref[cname],lvl[term],ref[term],gdensity[term]

#  if j==0 and sname == name[0]: # if we are at root, go for Venn diagram but this should be needed for last sample only
  if j==0: #alternative with obligatory Venn for all samples
  
   terms=sorted_terms

   ALL_GO=A|B|C
   print ' transcripts:', len(ALL_GO)
   plt.figure(figsize=(10,10))
   labels=(ref[terms[0]],ref[terms[1]],ref[terms[2]])
   v1=venn3([A,B,C],set_labels=labels, alpha =0.5)
   plt.title('Number of transcripts with GO:'+str(len(ALL_GO))+ '\nout of the total:'+str(len(r_list)))
   plt.show()
   plt.savefig(root_fname+'.venn_t_counts.pdf')
   plt.close()

   sets=[(A-B)-C, (B-A)-C, (A&B)-C, (C-A)-B, (A&C)-B, (B&C)-A, (A&B)&C]
   abut_list=[]
   for id_set in sets:
     abut_list.append(sum([abu_t[x][sname] for x in id_set]))
   print ' TPM:',int(sum(abut_list))
   #alt_abut=[abu_t[x][sname] for x in ALL_GO]
   #print 'all GO transcript abu:', sum(alt_abut)

   plt.figure(figsize=(10,10))
   v2=venn3(subsets=tuple(abut_list),set_labels=labels, alpha =0.5)
   plt.title('TPM covered by GO-assigned transcripts:'+ str(sum(abut_list)))
   plt.show()
   plt.savefig(root_fname+'.venn_t_abu.pdf')
   plt.close()

   # gene-level: convert sets of transcript IDs to sets of gene IDs
   # multi-gene transcripts will be counted in both sets
   # intersections will be relatively larger 
   # this should not make huge difference unless the assembly is really bad
   # if suspected, try: jaccard clip, pre-filtering reads on min coverage and re-assembly
   """
   AA=[]
   for tid in A:
    gid='_'.join(tid.split('_')[:-1])
    AA.append(gid)
   A=set(AA)
   
   BB=[]
   for tid in B:
    gid='_'.join(tid.split('_')[:-1])
    BB.append(gid)
   B=set(BB)
   
   CC=[]
   for tid in C:
    gid='_'.join(tid.split('_')[:-1])
    CC.append(gid)
   C=set(CC)
   
   """
   sets=[A,B,C]
   new_sets=[]
   for id_set in sets:
    mod_set=set()
    for each in id_set:
     g_id = '_'.join(each.split('_')[:-1])
     mod_set.add(g_id)
    new_sets.append(mod_set)
   A=new_sets[0]
   B=new_sets[1]
   C=new_sets[2]
   
   ALL_GO=A|B|C
   print ' genes:', len(ALL_GO)
   plt.figure(figsize=(10,10))
   v3=venn3([A,B,C], set_labels=labels, alpha =0.5)
   plt.title('Number of GO-assigned genes:'+ str(len(ALL_GO))+'\nout of the total ' +str(len(abu_d)))
   plt.show()
   plt.savefig(root_fname+'.venn_g_count.pdf')
   plt.close()

   sets=[(A-B)-C, (B-A)-C, (A&B)-C, (C-A)-B, (A&C)-B, (B&C)-A, (A&B)&C]
   abug_list=[]
   for id_set in sets:
    abug_list.append(sum([abu_d[x][sname] for x in id_set]))
   print ' TPM:', int(sum(abug_list))
   #alt_abug=[abu_d[x][sname] for x in ALL_GO]
   #print 'all GO gene abu:', sum(alt_abug)
   plt.figure(figsize=(10,10))
   v4=venn3(subsets=tuple(abug_list),set_labels=labels, alpha =0.5)
   plt.title('TPM covered by GO-assigned genes:'+ str(sum(abug_list)))
   plt.show()
   plt.savefig(root_fname+'.venn_g_abu.pdf')
   plt.close()
   
t11=time.time()
print t11-t10,'s for GO & Venn'
print'finished'

