import sys
import re
import os.path

file_path = sys.argv[1]
outname = file_path.split(".",1)[0]+'.gff'

with open(file_path, 'r') as f:
	lista_linii = f.readlines()
print 'reading', file_path
print 'appending annotations to', outname
ftype = 'rRNA'
source = "infernal"
phase = "."

with open(outname, 'a') as gff:
 for linia in lista_linii[2:-11]:
	kol = linia.split()
	start = kol[7]
	end = kol[8]
        score = kol[14]
        valid = kol[16] == '!'
	strand = kol[9]
        if strand == "-":
		start = kol[8]
		end = kol[7]
	aux = "product=" + kol[2]
	seqid = kol[0]
	gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + start + '\t' + end + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux +'\n'
        if valid:
          gff.write(gotowa_linia)

