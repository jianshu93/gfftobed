# gfftobed
Convert GFF3 to BED



This program takes an input genome annotation in GFF3 format and converts specific features to a 6 column BED format. It is capable of extracting any type of feature and any attribute in GFF3 format. It may also work with annotations in GTF. However, the 9th column of the GTF will be put in the 4th column of the BED file.

## Installation

Clone this repository, and `make` or `make clean`


## Usage


USAGE:
	gfftobed [options] <input_file_GFF3>

Extracts genomic coordinates of features from GFF3
-g/--gene			extract gene features in bed format
-e/--exons			extract exon features in bed format
-c/--cds			extract CDS features in bed format
-m/--mrna			extract mRNA features in bed format
-t/--tss			extract tss features in bed format
-f/--feature <feat>		extract custom features (i.e. ncRNA)

Extracts genomic coordinates of features with window around feature
-w/--window <int>		add <int> basepairs upstream and downstream of feature
-u/--upstream <int>		add <int> basepairs upstream/5' of feature
-d/--downstream <int>		add <int> basepairs downstream/3' of feature

-a/--attribute <string>		Specify attribute for name column ('note' by default)
-h/--help			Print this help message
  
  
