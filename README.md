# tritoconstrictor

Transcriptome annotation and filtering pipeline written in python and shell.</br>
Dependencies:</br>
-Python 2.7</br>
-Biopython</br>
-local copies of Pfam, Panther and UniRef90 databases</br>
-hmmer</br>
-infernal</br>
-MMseqs2</br>
-wise2</br>
-NCBI taxon database</br>
-InterPro mapping file</br>
-GO mapping file</br>
-Trinity assembler</br>
-hmm profiles of mitochondrial proteins from mitoconstrictor</br>
-cm profiles of mitochondrial rRNA from mitoconstrictor</br>
-cm profiles of rRNA subunits from Rfam</br>
</br>
Usage:</br>
-Prepare your assembly in fasta format </br>
-run database searches (see dbsearch.sh shell script for an example)</br>
-convert results to gff format (see the two prep_gff scripts)</br>
-prepare abundance estimate for each contig</br>
-run the main annotation script:</br>
  <i>python parse_gff3n.py </i> filename</br>
  where filename stands for your assembly</br>
  
