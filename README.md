# tritoconstrictor

Transcriptome annotation and filtering pipeline written in python and shell.<br>
Dependencies:<br>
-Python 2.7<br>
-Biopython<br>
-local copies of Pfam, Panther and UniRef90 databases<br>
-hmmer<br>
-infernal<br>
-MMseqs2<br>
-wise2<br>
-NCBI taxon database<br>
-InterPro mapping file<br>
-GO mapping file<br>
-Trinity assembler<br>
-hmm profiles of mitochondrial proteins from mitoconstrictor<br>
-cm profiles of mitochondrial rRNA from mitoconstrictor<br>
-cm profiles of rRNA subunits from Rfam<br>
<br>
Usage:<br>
-Prepare your assembly in fasta format<br>
-run database searches (see dbsearch.sh shell script for an example)<br>
-convert results to gff format (see the two prep_gff scripts)<br>
-prepare abundance estimate for each contig<br>
-run the main annotation script:<br>
  <i>python parse_gff3n.py </i> filename1 filename2<br>
  where filename1 stands for your assembly in fasta format and filename2 is the gff3 file<br>
 <br>
 Several files will be produced, including tbl and fsa files for GenBank submission, csv files with expression and GO data and gb file with annotated transcriptome.<br>
 It is expected that the first run will only give provisional results - further tuning is usually needed. This can be achieved by first editing bad taxon file, then re-running the annotation script. There is no need to repeat database searches. Once satisfied with the final filtering, re-run abundance estimate with the fsa file. Then the annotation script will produce more realistic abundance-based stats. Further tuning can be done if mitochondrial and rRNA loci are known - these can then be provided, substituting all mitochondrial and rRNA sequences.
  
