# gfftobed
Convert GFF3 to BED



This program takes an input genome annotation in GFF3 (1-based) format and converts specific features to a 6 column BED format (0-based). It is capable of extracting any type of feature and any attribute in GFF3 format. It may also work with annotations in GTF. However, the 9th column of the GTF will be put in the 4th column of the BED file.

## Installation

Clone this repository, and `make` or `make clean`


## Usage

USAGE:	`gfftobed [options] <input_file_GFF3>`

Option|Description
-------------------------|-------------------------------
-g/--gene|extract gene features in bed format
-e/--exons|extract exon features in bed format
-c/--cds|extract CDS features in bed format
-m/--mrna|extract mRNA features in bed format
-t/--tss|extract tss features in bed format
-f/--feature `<feat>`|extract custom features (i.e. ncRNA)
-w/--window `<int>`|add `<int>` basepairs upstream and downstream of feature
-u/--upstream `<int>`|add `<int>` basepairs upstream/5' of feature
-d/--downstream `<int>`|add `<int>` basepairs downstream/3' of feature
-a/--attribute `<string>`|Specify attribute for name column ('note' by default)
-h/--help|Print help message


## Troubleshooting

Check the validity of the input GFF3. If there is a formatting error in a particular line, it will either be skipped or will cause a problem. If there is a problem with the 4th column of the BED file (Name), try changing the option argument passed to -a/--attribute. If the input GFF3 file contains "note=" and you do not want the note attribute and the others are not working, try entering a nonsense string of characters to the option argument for -a/--attribute (example: -a lkalfkdjlk). This version was designed to work with a specific GFF3, future versions will be more robust.
