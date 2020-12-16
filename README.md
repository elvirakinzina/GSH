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
Pipeline for identification for novel human Genomic Safe Harbor (GSH) sites.
The following criteria were used to computationally predict novel GSH:
- 50kb away from known genes
- 300kb away from known oncogenes
- 300kb away from microRNAs
- Outside lncRNAs, tRNAs, centromeres, telomeres, genomic gaps

## Features
Identified GSH sites allow for stable and safe expression of genes of interest. The universality of the GSH selection criteria allow for the use of these sites in any cell type.

## Prerequisites
- wig2bed and gtf2bed from BEDOPS https://bedops.readthedocs.io/en/latest/content/installation.html#installation
- bedtools https://bedtools.readthedocs.io/en/latest/content/installation.html

## Data
- GENCODE annotation of the human genome version 24: comprehensive gene annotation, long non-coding RNA gene annotation, and predicted tRNA genes
https://www.gencodegenes.org/human/ 
- Cancer Gene Census (CGC) genes from tier 1 (extensive evidence of association with cancer available) and tier 2 (strong indications of the association exist) https://cancer.sanger.ac.uk/census
- microRNAs from MirGeneDB (based on Fromm et al, A Uniform System For The Annotation Of Human microRNA Genes And The Evolution Of The Human microRNAome. Annu Rev Genet. 2015, 49: 213–242) https://mirgenedb.org
- Genome sequence, GRCh38.p5 assembly
- Centromeric regions http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
- Gaps in the genome sequence which include telomeres http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz

## Usage

Usage:
```bash
  ./predict_gsh.sh [-genes] [-oncogenes] [-micrornas] [-trnas] [-lncrnas] [-enhancers] [-enhancers] [-centromeres] [] [-dist_from_genes] [-dist_from_oncogenes [-dist_from_micrornas] [-dist_from_micrornas] [-dist_from_trnas] [-dist_from_lncrnas] [-dist_from_enhancers] [-dist_from_centromeres] [-dist_from_gaps] [-h|--help]	
```
  
Options:
```bash
	-genes: Whether to exclude regions with and around genes (default=true)
	-oncogenes: Whether to exclude regions with and around oncogenes (default=true)
	-micrornas: Whether to exclude regions with and around microRNAs (default=true)
	-trnas: Whether to exclude regions with and around tRNAs (default=true)
	-lncrnas: Whether to exclude regions with and around lncRNAs (default=true)
	-enhancers: Whether to exclude regions with and around enhancers (default=false)
	-centromeres: Whether to exclude regions with and around centromeres (default=true)
	-gaps: Whether to exclude regions with and around gaps (default=true)
	-dist_from_genes: Minimal distance from any safe harbor to any gene in bp (default=50000)
	-dist_from_oncogenes: Minimal distance from any safe harbor to any oncogene in bp (default=300000)
	-dist_from_micrornas: Minimal distance from any safe harbor to any microRNA in bp (default=300000)
	-dist_from_trnas: Minimal distance from any safe harbor to any tRNA in bp (default=0)
	-dist_from_lncrnas: Minimal distance from any safe harbor to any long-non-coding RNA in bp (default=0)
	-dist_from_enhancers: Minimal distance from any safe harbor to any enhancer in bp (default=300000)
	-dist_from_centromeres: Minimal distance from any safe harbor to any centromere in bp (default=0)
	-dist_from_gaps: Minimal distance from any safe harbor to any gaps in bp (default=0)
	-h, --help: Prints help
 ```
 
Example:
```bash
./predict_gsh.sh -dist_from_trnas 0 -dist_from_gaps 0 -dist_from_centromeres 0 -enhancers false -dist_from_lncrnas 0
Distance from genes = 50000 bp
Distance from microRNAs = 300000 bp
Distance from tRNAs = 0 bp
Distance from lncRNAs = 0 bp
Distance from oncogenes = 300000 bp
Distance from centromeres = 0 bp
Distance from gaps = 0 bp
Merging all genomic regions to avoid
Obtaining genomic coordinates and sequences of safe harbors
 ```
 
## Reference
Aznauryan et al. (2020), Discovery and validation of novel human genomic safe harbor sites for gene and cell therapies. Cell Genomics. In review

## License

