#!/bin/bash

###################################################################################
#
#	        FILE: predict_gsh.sh
#
#	       USAGE: ./predict_gsh.sh
#
#	 DESCRIPTION: Find all universal genomic safe harbors in the human genome     
#
#	      AUTHOR: Elvira Kinzina, elvira@mit.edu
#
#	     CREATED: 12/10/2020 10:02 AM              
#              
###################################################################################

#=======================================================================
#
#							 Arguments
#
#=======================================================================

# THE DEFAULTS INITIALIZATION
# ./predict_gsh.sh -dist_from_trnas 0 -dist_from_gaps 0 -dist_from_centromeres 0 -enhancers false -dist_from_lncrnas 0
genes=true
oncogenes=true
micrornas=true
t_cell_specific_micrornas=true
trnas=true
lncrnas=true
enhancers=true
centromeres=true
gaps=true
dist_from_genes=50000
dist_from_oncogenes=300000
dist_from_micrornas=300000
dist_from_trnas=150000
dist_from_lncrnas=150000
dist_from_enhancers=20000
dist_from_centromeres=150000
dist_from_gaps=150000

print_help()
{
	printf '%s\n' "Argbash is an argument parser generator for Bash."
	printf 'Usage: %s [-dist_from_genes] [-dist_from_oncogenes [-dist_from_micrornas] [-dist_from_micrornas] [-dist_from_trnas] [-dist_from_lncrnas] [-dist_from_enhancers] [-dist_from_centromeres] [-dist_from_gaps] [-h|--help]'
	printf '\t%s\n' "-genes: Whether to exclude regions with and around genes (default=true)"
	printf '\t%s\n' "-oncogenes: Whether to exclude regions with and around oncogenes (default=true)"
	printf '\t%s\n' "-micrornas: Whether to exclude regions with and around microRNAs (default=true)"
	printf '\t%s\n' "-trnas: Whether to exclude regions with and around tRNAs (default=true)"
	printf '\t%s\n' "-lncrnas: Whether to exclude regions with and around lncRNAs (default=true)"
	printf '\t%s\n' "-enhancers: Whether to exclude regions with and around enhancers (default=true)"
	printf '\t%s\n' "-centromeres: Whether to exclude regions with and around centromeres (default=true)"
	printf '\t%s\n' "-gaps: Whether to exclude regions with and around gaps (default=true)"
	printf '\t%s\n' "-dist_from_genes: Minimal distance from any safe harbor to any gene in bp (default=50000)"
	printf '\t%s\n' "-dist_from_oncogenes: Minimal distance from any safe harbor to any oncogene in bp (default=300000)"
	printf '\t%s\n' "-dist_from_micrornas: Minimal distance from any safe harbor to any microRNA in bp (default=300000)"
	printf '\t%s\n' "-dist_from_trnas: Minimal distance from any safe harbor to any tRNA in bp (default=150000)"
	printf '\t%s\n' "-dist_from_lncrnas: Minimal distance from any safe harbor to any long-non-coding RNA in bp (default=150000)"
	printf '\t%s\n' "-dist_from_enhancers: Minimal distance from any safe harbor to any enhancer in bp (default=20000)"
	printf '\t%s\n' "-dist_from_centromeres: Minimal distance from any safe harbor to any centromere in bp (default=150000)"
	printf '\t%s\n' "-dist_from_gaps: Minimal distance from any safe harbor to any gaps in bp (default=150000)"
	printf '\t%s\n' "-h, --help: Prints help"
}


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    print_help
    exit 0
    ;;
    -genes)
    genes="$2"
    shift # past argument
    shift # past value
    ;;
    -oncogenes)
    oncogenes="$2"
    shift # past argument
    shift # past value
    ;;
    -micrornas)
    micrornas="$2"
    shift # past argument
    shift # past value
    ;;
    -t_cell_specific_micrornas)
    t_cell_specific_micrornas="$2"
    shift # past argument
    shift # past value
    ;;
    -trnas)
    trnas="$2"
    shift # past argument
    shift # past value
    ;;
    -lncrnas)
    lncrnas="$2"
    shift # past argument
    shift # past value
    ;;
    -enhancers)
    enhancers="$2"
    shift # past argument
    shift # past value
    ;;
    -centromeres)
    centromeres="$2"
    shift # past argument
    shift # past value
    ;;
    -gaps)
    gaps="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_genes)
    dist_from_genes="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_oncogenes)
    dist_from_oncogenes="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_micrornas)
    dist_from_micrornas="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_trnas)
    dist_from_trnas="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_lncrnas)
    dist_from_lncrnas="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_enhancers)
    dist_from_enhancers="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_centromeres)
    dist_from_centromeres="$2"
    shift # past argument
    shift # past value
    ;;
    -dist_from_gaps)
    dist_from_gaps="$2"
    shift # past argument
    shift # past value
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#=======================================================================

DIRECTORY=tmp
if [ -d "$DIRECTORY" ]; then
  rm -r $DIRECTORY
fi

mkdir tmp	# Create temporary folder for intermediate files

#=======================================================================
#
#						     Genes
#
#=======================================================================

if [ "$genes" = true ] ; then
	echo "Distance from genes = ${dist_from_genes}" bp

	mkdir tmp/genes
	dir=tmp/genes

	# get gene annotation from GENCODE
#	less gencode.v24.annotation.gtf | grep "\tgene\t" >> ${dir}/gencode_gene_anotation.gtf

#	awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' ${dir}/gencode_gene_anotation.gtf >> ${dir}/gencode_v24_annotation_genes_transcript_id.gtf

	# get genomic coordinates of all annotated genes in the human genome in the BED format
#	gtf2bed < ${dir}/gencode_v24_annotation_genes_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_annotation_genes.bed

	# get genomic regions of length ${dist_from_genes} base pairs flanking genes from both sides
	bedtools slop -b ${dist_from_genes} -i data/gencode_v24_annotation_genes.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_annotation_genes_with_flanks.bed

	# merge regions containing genes and their flanking regions
	sortBed -i ${dir}/gencode_v24_annotation_genes_with_flanks.bed >> ${dir}/gencode_v24_annotation_genes_with_flanks_sorted.bed

	bedtools merge -i ${dir}/gencode_v24_annotation_genes_with_flanks_sorted.bed >> ${dir}/gencode_v24_annotation_genes_with_flanks_merged.bed
fi

#=======================================================================
#
#							 Oncogenes
#
#=======================================================================

if [ "$oncogenes" = true ] ; then
	echo "Distance from oncogenes = ${dist_from_oncogenes}" bp

	mkdir tmp/oncogenes
	dir=tmp/oncogenes

	# get GENCODE gene annotation for oncogenes from COSMIC (Cancer Gene Census)
#	grep -f data/Oncogenes_id_list.txt tmp/genes/gencode_v24_annotation_genes_transcript_id.gtf >> tmp/oncogenes/gencode_oncogenes_annotation_transcript_id.gtf

	# get genomic coordinates of the oncogenes in the BED format
#	gtf2bed < ${dir}/gencode_oncogenes_annotation_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_annotation_oncogenes.bed

	# get genomic regions of length ${dist_from_oncogenes} base pairs flanking oncogenes from both sides
	bedtools slop -b ${dist_from_oncogenes} -i data/gencode_v24_annotation_oncogenes.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_annotation_oncogenes_with_flanks.bed

	# merge regions containing oncogenes and their flanking regions
	sortBed -i ${dir}/gencode_v24_annotation_oncogenes_with_flanks.bed >> ${dir}/gencode_v24_annotation_oncogenes_with_flanks_sorted.bed

	bedtools merge -i ${dir}/gencode_v24_annotation_oncogenes_with_flanks_sorted.bed >> ${dir}/gencode_v24_annotation_oncogenes_with_flanks_merged.bed
fi

#=======================================================================
#
#						   MicroRNAs
#
#=======================================================================

if [ "$micrornas" = true ] ; then
	echo "Distance from microRNAs = ${dist_from_micrornas}" bp
	mkdir tmp/micrornas
	dir=tmp/micrornas
	
	if [ "$t_cell_specific_micrornas" = true ] ; then
		# convert wig files from ENCODE to BED files
		wig2bed --keep-header < data/ENCFF417SFH.wig > ${dir}/ENCFF417SFH_bedops.bed
		wig2bed --keep-header < data/ENCFF850OSS.wig > ${dir}/ENCFF850OSS_bedops.bed

		# remove lines containing '_header'
		grep -v '_header' ${dir}/ENCFF417SFH_bedops.bed >> ${dir}/ENCFF417SFH.bed
		grep -v '_header' ${dir}/ENCFF850OSS_bedops.bed >> ${dir}/ENCFF850OSS.bed

		# exclude regions aligned to the EBV genome (chrEBV) 
		less ${dir}/ENCFF417SFH.bed | grep -v chrEBV >> ${dir}/ENCFF417SFH_no_chrEBV.bed 
		less ${dir}/ENCFF850OSS.bed | grep -v chrEBV >> ${dir}/ENCFF850OSS_no_chrEBV.bed

		# remove lines containing '_header'
		grep -v '_header' ${dir}/ENCFF417SFH_no_chrEBV.bed >> ${dir}/ENCFF417SFH_no_chrEBV2.bed
		grep -v '_header' ${dir}/ENCFF850OSS_no_chrEBV.bed >> ${dir}/ENCFF850OSS_no_chrEBV2.bed

		# get genomic regions of length ${dist_from_micrornas} base pairs flanking microRNAs from both sides
		bedtools flank -b ${dist_from_micrornas} -i ${dir}/ENCFF417SFH_no_chrEBV2.bed -g data/chromInfo_hg38.txt >> ${dir}/ENCFF417SFH_flanks.bed
		bedtools flank -b ${dist_from_micrornas} -i ${dir}/ENCFF850OSS_no_chrEBV2.bed -g data/chromInfo_hg38.txt >> ${dir}/ENCFF850OSS_flanks.bed

		# merge regions containing microRNAs and their flanking regions
		cat ${dir}/ENCFF850OSS_no_chrEBV.bed ${dir}/ENCFF850OSS_flanks.bed >> ${dir}/ENCFF850OSS_with_flanks.bed
		cat ${dir}/ENCFF417SFH_no_chrEBV.bed ${dir}/ENCFF417SFH_flanks.bed >> ${dir}/ENCFF417SFH_with_flanks.bed

		sortBed -i ${dir}/ENCFF417SFH_with_flanks.bed >> ${dir}/ENCFF417SFH_with_flanks_sorted.bed
		sortBed -i ${dir}/ENCFF850OSS_with_flanks.bed >> ${dir}/ENCFF850OSS_with_flanks_sorted.bed

		bedtools merge -i ${dir}/ENCFF417SFH_with_flanks_sorted.bed >> ${dir}/ENCFF417SFH_with_flanks_merged.bed
		bedtools merge -i ${dir}/ENCFF850OSS_with_flanks_sorted.bed >> ${dir}/ENCFF850OSS_with_flanks_merged.bed
	fi
	
	if [ "$t_cell_specific_micrornas" = false ] ; then
		# get genomic regions of length ${dist_from_micrornas} base pairs flanking microRNAs from both sides
		bedtools slop -b ${dist_from_micrornas} -i data/hsa-all.bed -g data/chromInfo_hg38.txt >> ${dir}/Micrornas_with_flanks.bed

		# merge regions containing microRNAs and their flanking regions
		sortBed -i ${dir}/Micrornas_with_flanks.bed >> ${dir}/Micrornas_with_flanks_sorted.bed

		bedtools merge -i ${dir}/Micrornas_with_flanks_sorted.bed >> ${dir}/Micrornas_with_flanks_merged.bed
	fi
	
fi

#=======================================================================
#
#					     Long-non-coding RNAs
#
#=======================================================================

if [ "$lncrnas" = true ] ; then
	echo "Distance from lncRNAs = ${dist_from_lncrnas}" bp

	mkdir tmp/lncrnas
	dir=tmp/lncrnas

	# get lncRNA annotation from GENCODE
#	awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v24.long_noncoding_RNAs.gtf >> ${dir}/gencode_v24_long_noncoding_RNAs_transcript_id.gtf

#	gtf2bed < ${dir}/gencode_v24_long_noncoding_RNAs_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_long_noncoding_RNAs.bed

	# get genomic regions of length ${dist_from_lncrnas} base pairs flanking lncRNAs from both sides
	bedtools slop -b ${dist_from_lncrnas} -i data/gencode_v24_long_noncoding_RNAs.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks.bed

	# merge regions containing lncRNAs and their flanking regions
	sortBed -i ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks.bed >> ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_sorted.bed

	bedtools merge -i ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_sorted.bed >> ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_merged.bed
fi

#=======================================================================
#
#						     tRNAs
#
#=======================================================================

if [ "$trnas" = true ] ; then
	echo "Distance from tRNAs = ${dist_from_trnas}" bp

	mkdir tmp/trnas
	dir=tmp/trnas

	# get tRNA annotation from GENCODE
	gtf2bed < gencode.v24.tRNAs.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_tRNAs.bed

	# get genomic regions of length ${dist_from_trnas} base pairs flanking tRNAs from both sides
	bedtools slop -b ${dist_from_trnas} -i ${dir}/gencode_v24_tRNAs.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_tRNAs_with_flanks.bed

	# merge regions containing tRNAs and their flanking regions
	sortBed -i ${dir}/gencode_v24_tRNAs_with_flanks.bed >> ${dir}/gencode_v24_tRNAs_with_flanks_sorted.bed

	bedtools merge -i ${dir}/gencode_v24_tRNAs_with_flanks_sorted.bed >> ${dir}/gencode_v24_tRNAs_with_flanks_merged.bed
	
fi

#=======================================================================
#
#							 Enhancers
#
#=======================================================================

if [ "$enhancers" = true ] ; then
	echo "Distance from enhancers = ${dist_from_enhancers}" bp

	mkdir tmp/enhancers
	dir=tmp/enhancers
	
	# get genomic regions of length ${dist_from_enhancers} base pairs flanking enhancers from both sides
	bedtools slop -b ${dist_from_enhancers} -i data/All_human_enhancers_sorted.bed -g data/chromInfo_hg38.txt >> ${dir}/All_human_enhancers_with_flanks.bed

	# merge regions containing enhancers and their flanking regions
	sortBed -i ${dir}/All_human_enhancers_with_flanks.bed >> ${dir}/All_human_enhancers_with_flanks_sorted.bed

	bedtools merge -i ${dir}/All_human_enhancers_with_flanks_sorted.bed >> ${dir}/All_human_enhancers_with_flank_merged.bed

fi

#=======================================================================
#
#    						Centromeres
#
#=======================================================================

if [ "$centromeres" = true ] ; then
	echo "Distance from centromeres = ${dist_from_centromeres}" bp

	mkdir tmp/centromeres
	dir=tmp/centromeres

	# get genomic coordinates of centromeres in the BED format
	less data/hgTables_centromeric_regions_38.txt | tail -n +3 >> ${dir}/hgTables_centromeric_regions_38.bed
	less ${dir}/hgTables_centromeric_regions_38.bed | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/Centromeres.bed 

	# get genomic regions of length ${dist_from_centromeres} base pairs flanking centromeres from both sides
	bedtools slop -b ${dist_from_centromeres} -i ${dir}/Centromeres.bed -g data/chromInfo_hg38.txt >> ${dir}/Centromeres_with_flanks.bed

	# merge regions containing centromeres and their flanking regions
	sortBed -i ${dir}/Centromeres_with_flanks.bed >> ${dir}/Centromeres_with_flanks_sorted.bed

	bedtools merge -i ${dir}/Centromeres_with_flanks_sorted.bed >> ${dir}/Centromeres_with_flanks_merged.bed
fi

#=======================================================================
#
#   	Gaps: telomeres and other regions that cannot be sequenced
#
#=======================================================================

if [ "$gaps" = true ] ; then
	echo "Distance from gaps = ${dist_from_gaps}" bp

	mkdir tmp/gaps
	dir=tmp/gaps

	# get genomic coordinates of gaps in the BED format
	less data/hgTables_gaps.txt | tail -n +2 | cut -f 2,3,4 >> ${dir}/hgTables_gaps.bed

	# get genomic regions of length ${dist_from_centromeres} base pairs flanking gaps from both sides
	bedtools slop -b ${dist_from_gaps} -i ${dir}/hgTables_gaps.bed -g data/chromInfo_hg38.txt >> ${dir}/Gaps_with_flanks.bed

	# merge regions containing gaps and their flanking regions
	sortBed -i ${dir}/Gaps_with_flanks.bed >> ${dir}/Gaps_with_flanks_sorted.bed

	bedtools merge -i ${dir}/Gaps_with_flanks_sorted.bed >> ${dir}/Gaps_with_flanks_merged.bed
fi

#=======================================================================
#
#		 		  Union of all genomic regions to avoid
#
#=======================================================================

echo "Merging all genomic regions to avoid"

mkdir tmp/merge
dir=tmp/merge

# intersect of regions to avoid together
if [ "$genes" = true ] ; then
	cat tmp/genes/gencode_v24_annotation_genes_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$genes" = true ] ; then
	cat tmp/oncogenes/gencode_v24_annotation_oncogenes_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$micrornas" = true ] ; then
	if [ "$t_cell_specific_micrornas" = true ] ; then
		cat tmp/micrornas/ENCFF850OSS_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
	fi
	
	if [ "$t_cell_specific_micrornas" = false ] ; then
		cat tmp/micrornas/Micrornas_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
	fi
fi

if [ "$trnas" = true ] ; then
	cat tmp/trnas/gencode_v24_tRNAs_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$lncrnas" = true ] ; then
	cat tmp/lncrnas/gencode_v24_long_noncoding_RNAs_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$enhancers" = true ] ; then
	cat tmp/enhancers/All_human_enhancers_with_flank_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$centromeres" = true ] ; then
	cat tmp/centromeres/Centromeres_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi

if [ "$gaps" = true ] ; then
	cat tmp/gaps/Gaps_with_flanks_merged.bed >> ${dir}/Regions_to_avoid.bed
fi


sortBed -i ${dir}/Regions_to_avoid.bed >> ${dir}/Regions_to_avoid_sorted.bed

bedtools merge -i ${dir}/Regions_to_avoid_sorted.bed >> ${dir}/Regions_to_avoid_merged.bed

#=======================================================================
#
#							Safe harbors
#
#=======================================================================

echo "Obtaining genomic coordinates and sequences of safe harbors"

mkdir tmp/safe_harbors
dir=tmp/safe_harbors

# substract all regions to avoid from the whole genome
bedtools subtract -a data/all_chrom.bed -b tmp/merge/Regions_to_avoid_merged.bed >> ${dir}/Safe_harbors_with_alt.bed

# exclude pseudo-chromosomes and alterative loci
grep -v '_' ${dir}/Safe_harbors_with_alt.bed >> ${dir}/Safe_harbors.bed

FILE=Safe_harbors.bed
if [[ -f "$FILE" ]]; then
    rm $FILE
fi

sortBed -i ${dir}/Safe_harbors.bed >> Safe_harbors.bed

FILE=Safe_harbors.fasta
if [[ -f "$FILE" ]]; then
    rm $FILE
fi

# get sequences of those regions
bedtools getfasta -fi GRCh38.p5.genome.fa -bed Safe_harbors.bed >> Safe_harbors.fasta
