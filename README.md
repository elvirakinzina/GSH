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
- Enhancers from ROADMAP https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz



## Usage

Running:
```bash
Usage:  ./predict_gsh.sh [-genes] [-oncogenes] [-micrornas] [-trnas] [-lncrnas] [-enhancers] [-enhancers] [-centromeres] [] [-dist_from_genes] [-dist_from_oncogenes [-dist_from_micrornas] [-dist_from_micrornas] [-dist_from_trnas] [-dist_from_lncrnas] [-dist_from_enhancers] [-dist_from_centromeres] [-dist_from_gaps] [-h|--help]	
```
  
Options:
```bash
	-genes: Whether to exclude regions with and around genes (default=true)
	-oncogenes: Whether to exclude regions with and around oncogenes (default=true)
	-micrornas: Whether to exclude regions with and around microRNAs (default=true)
	-trnas: Whether to exclude regions with and around tRNAs (default=true)
	-lncrnas: Whether to exclude regions with and around lncRNAs (default=true)
	-enhancers: Whether to exclude regions with and around enhancers (default=true)
	-centromeres: Whether to exclude regions with and around centromeres (default=true)
	-gaps: Whether to exclude regions with and around gaps (default=true)
	-dist_from_genes: Minimal distance from any safe harbor to any gene in bp (default=50000)
	-dist_from_oncogenes: Minimal distance from any safe harbor to any oncogene in bp (default=300000)
	-dist_from_micrornas: Minimal distance from any safe harbor to any microRNA in bp (default=300000)
	-dist_from_trnas: Minimal distance from any safe harbor to any tRNA in bp (default=300000)
	-dist_from_lncrnas: Minimal distance from any safe harbor to any long-non-coding RNA in bp (default=300000)
	-dist_from_enhancers: Minimal distance from any safe harbor to any enhancer in bp (default=300000)
	-dist_from_centromeres: Minimal distance from any safe harbor to any centromere in bp (default=50000)
	-dist_from_gaps: Minimal distance from any safe harbor to any gaps in bp (default=50000)
	-h, --help: Prints help
 ```
 
## Reference
Aznauryan et al. (2020), Bioinformatic screening and phenotypic validation of novel genomic safe harbor sites by targeted integration and transcriptome profiling

## License

