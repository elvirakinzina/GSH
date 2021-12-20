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
- 300kb away from microRNAs, centromeres, telomeres, genomic gaps
- 150kb away from lncRNAs, tRNAs
- 20kb away from enhancers

## Features
Produced genomic coordinates serve as computationally predicted GSH sites based on previously established as well as newly introduced criteria. All criteria used are universal for all cell types.

## Prerequisites
- wig2bed and gtf2bed from BEDOPS https://bedops.readthedocs.io/en/latest/content/installation.html#installation
- bedtools https://bedtools.readthedocs.io/en/latest/content/installation.html

## Data
- GENCODE annotation of the human genome version 24: comprehensive gene annotation, long non-coding RNA gene annotation, and predicted tRNA genes
https://www.gencodegenes.org/human/ 
- Cancer Gene Census (CGC) genes from tier 1 (extensive evidence of association with cancer available) and tier 2 (strong indications of the association exist) https://cancer.sanger.ac.uk/census
- microRNAs from MirGeneDB (based on Fromm et al, A Uniform System For The Annotation Of Human microRNA Genes And The Evolution Of The Human microRNAome. Annu Rev Genet. 2015, 49: 213–242) https://mirgenedb.org
- Enhancers were obtained from EnhacerAtlas 2.0 http://www.enhanceratlas.org/downloadv2.php. Coordinates for human enhancers are given for the hg19 genome assembly. To convert to the GRCh38/hg38 version, use the LiftOver tool from the UCSC genome browser https://genome.ucsc.edu/cgi-bin/hgLiftOver
- Centromeric regions http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
- Gaps in the genome sequence which include telomeres http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
- Genome sequence, GRCh38.p5 assembly https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_5fe9f269a15a7a0665583144&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true

## Usage

Usage:
```bash
  ./predict_gsh.sh [-genes] [-oncogenes] [-micrornas] [-trnas] [-lncrnas] [-enhancers] [-centromeres] [-gaps] [-dist_from_genes] [-dist_from_oncogenes] [-dist_from_micrornas] [-dist_from_trnas] [-dist_from_lncrnas] [-dist_from_enhancers] [-dist_from_centromeres] [-dist_from_gaps] [-h|--help]	
```
  
Options:
```bash
	-genes: Whether to exclude regions with and around genes (default=true)
	-oncogenes: Whether to exclude regions with and around oncogenes (default=true)
	-micrornas: Whether to exclude regions with and around microRNAs (default=true)
	-trnas: Whether to exclude regions with and around tRNAs (default=true)
	-lncrnas: Whether to exclude regions with and around lncRNAs (default=true)
	-centromeres: Whether to exclude regions with and around centromeres (default=true)
	-gaps: Whether to exclude regions with and around gaps (default=true)
	-enhancers: Whether to exclude enhancer regions (default=true)

	-dist_from_genes: Minimal distance from any safe harbor to any gene in bp (default=50000)
	-dist_from_oncogenes: Minimal distance from any safe harbor to any oncogene in bp (default=300000)
	-dist_from_micrornas: Minimal distance from any safe harbor to any microRNA in bp (default=300000)
	-dist_from_trnas: Minimal distance from any safe harbor to any tRNA in bp (default=150000)
	-dist_from_lncrnas: Minimal distance from any safe harbor to any long-non-coding RNA in bp (default=150000)
	-dist_from_enhancers: Minimal distance from any safe harbor to any enhancer in bp (default=20000)
	-dist_from_centromeres: Minimal distance from any safe harbor to any centromere in bp (default=300000)
	-dist_from_gaps: Minimal distance from any safe harbor to any gaps in bp (default=300000)
	-h, --help: Prints help
 ```
 
Running with the default parameters:
```bash
./predict_gsh.sh ```
Distance from genes = 50000 bp
Distance from oncogenes = 300000 bp
Distance from microRNAs = 300000 bp
Distance from tRNAs = 150000 bp
Distance from lncRNAs = 150000 bp
Distance from enhancers = 20000 bp
Distance from centromeres = 300000 bp
Distance from gaps = 300000 bp
Merging all genomic regions to avoid
Obtaining genomic coordinates and sequences of safe harbors
 
Running with modified parameters:
```bash
./predict_gsh.sh -trnas false -dist_from_lncrnas 50000 -dist_from_enhancers 0 -dist_from_centromeres 0 -dist_from_gaps 0
```
Distance from genes = 50000 bp
Distance from oncogenes = 300000 bp
Distance from microRNAs = 300000 bp
Distance from lncRNAs = 50000 bp
Distance from enhancers = 0 bp
Distance from centromeres = 0 bp
Distance from gaps = 0 bp
Merging all genomic regions to avoid
Obtaining genomic coordinates and sequences of safe harbors
 
The output is two files: Safe_harbors.bed that has genomic coordinates of all regions potentially containing safe harbors and Safe_harbors.fasta contains sequences of these regions.

## Reference
Aznauryan et al. (2020), Discovery and validation of novel human genomic safe harbor sites for gene and cell therapies. Cell Genomics. In review

