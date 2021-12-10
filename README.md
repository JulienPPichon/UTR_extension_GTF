This script has been developped in the frame of the M1 internship. 
It extends the transcripts' 3'UTRs of GTF file to have a better mapping of 
library with a high number of reads in 3'UTR. (e.g. ScRNA-seq analysis).

usage: gene_extension.py [-h] [-i FILENAME.gtf] [-o OUTPUT] [-l EXTENT_LENGTH]
                         [-d GENE_DISTANCE] [-e] [-c CHROMOSOME_SIZE]

optional arguments:

  -h, --help            show this help message and exit
  
  -i FILENAME, --filename FILENAME
                        take filename.gtf as input
                        
  -o OUTPUT, --output OUTPUT
                        name of the gtf output file with extended genes
                        
  -l EXTENT_LENGTH, --extent_length EXTENT_LENGTH
                        Length of last exon transcript extension
                        
  -d GENE_DISTANCE, --gene_distance GENE_DISTANCE
                        minimum distance separing 2 transcripts of different
                        genes
                        
  -e, --exon_number     add exon number at the end of the attributes
  
  -c CHROMOSOME_SIZE, --chromosome_size CHROMOSOME_SIZE
                        Use chromosome size file to extend last gene
