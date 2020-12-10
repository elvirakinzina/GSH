# GSH - Genomic Safe Harbors

## Index

* [Description](#description)
* [Features](#features)
* [Prerequisites](#prerequisites)
* [Data](#usage)
* [Usage](#usage)
* [Reference](#reference)
* [License](#license)

## Description

## Features

## Prerequisites
- bigwig
- wig2bed
- bedtools

## Data
- microRNAs from Encode (microRNA-seq of CD4-positive, alpha-beta memory T cell, ENCSR199AQA, plus and minus strand)

## Usage

Commands used:

# bedtools sorted
```bash
$ bedtools intersect \
           -a ccds.exons.bed -b aln.bam.bed \
           -c \
           -sorted
```

Usage:

```bash
  bedtools (-h | --help)
```
  
Options:
```bash
  -h --help     Show this screen.
  --version     Show version.
 ```
 
## Reference
Aznauryan et al. (2020), Bioinformatic screening and phenotypic validation of novel genomic safe harbor sites by targeted integration and transcriptome profiling

## License

