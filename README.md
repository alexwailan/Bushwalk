# Bushwalk

## Introduction
Lodestone provides variant calls as vcf files. Bushwalk takes the Lodestone output and prepares the data of each sample for input into the snippy-core module of Snippy.  We're going on a bushwalk!

## Dependencies
------

-   samtools
-   bcftools


## Inputs
------

-   Lodestone outputs vcf files in their own 

## Usage

usage: potplant.py tree [options]

Run PotPlant.py A program to root and collpase branches of zero length.

positional arguments:
  FILE                  The tree file in newick format.

optional arguments:
  -h, --help            show this help message and exit
  -r FILE [FILE ...], --root FILE [FILE ...]
                        The isolate for rooting the tree.
  -p PREFIX [PREFIX ...], --prefix PREFIX [PREFIX ...]
                        The prefix for the output.
  -d [N], --dirpath [N]
                        Input directory containing the tree. End with a
                        forward slash. Eg. /temp/fasta/
  -o [N], --outdir [N]  Output directory. End with a forward slash. Eg.
                        /temp/fasta/; Default to use current directory.
----------
