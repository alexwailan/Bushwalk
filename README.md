# Bushwalk

## Introduction
Sometimes you have variant calls as vcf files. Bushwalk takes the vcf files and prepares the data of each sample for input into the snippy-core module of Snippy with the intended use to creat a multisequence alignment.  We're going on a bushwalk!

## Dependencies
------

-   samtools
-   bcftools


## Inputs
------

-   Lodestone outputs vcf files in their own directory under the file structure ID/lodestone/
-   Reference file paired with the vcf file: Needs to be in FASTA format and be only the chromosome
-   A file containing all the IDs of the samples corresponding. Each sample ID will be on a new line.

## Usage
```

usage: bushwalk.py [options] reference ids

Run Bushwalk. A program to parse the output of Lodestone for input into
snippy-core

positional arguments:
  N                     Provide the reference for snippy [Required]
  N                     IDs of paired reads files [Required]

optional arguments:
  -h, --help            show this help message and exit
  -d [N], --dirpath [N]
                        Input directory of Lodestone results. End with a
                        forward slash. Eg. /temp/fasta/ [Required]
  -o [N], --outdir [N]  Output directory. End with a forward slash. Eg.
                        /temp/fasta/; Default to use current directory.

```
## Output
Bushwalk will generate:
- A new VCF file containing only SNV calls
- A SNV consensus fasta file
