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
Finding all genomic safe harbors in the human genome that are universal for most cell types   

## Features

## Prerequisites
- wig2bed and gtf2bed from BEDOPS https://bedops.readthedocs.io/en/latest/content/installation.html#installation
- bedtools https://bedtools.readthedocs.io/en/latest/content/installation.html

## Data
- GENCODE annotation of the human genome version 24
https://www.gencodegenes.org/human/ 
- microRNAs from ENCODE (microRNA-seq of CD4-positive, alpha-beta memory T cell, ENCSR199AQA, plus and minus strand)
- t-RNAs


## Usage

```bash
$ ./predict_gsh.sh -dist_from_genes 50000 -dist_from_oncogenes 300000 -dist_from_micrornas 300000 -dist_from_trnas 300000 -dist_from_lncrnas 300000

```
  
Options:
```bash
  -h --help     Show this screen.
  --version     Show version.
	-dist_from_genes \
	-dist_from_oncogenes \
	-dist_from_micrornas \
	-dist_from_trnas \
	-dist_from_lncrnas
 ```
 
## Reference
Aznauryan et al. (2020), Bioinformatic screening and phenotypic validation of novel genomic safe harbor sites by targeted integration and transcriptome profiling

## License

