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
#
#				    Arguments
#
#==================================================================================

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
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
    -dist_from_enhancers)
    dist_from_enhancers="$2"
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
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#==================================================================================

mkdir tmp	# Create temporary folder for intermediate files

#==================================================================================
#
#				     Genes
#
#==================================================================================

echo "Distance from genes = ${dist_from_genes}" bp

mkdir tmp/genes
dir=tmp/genes

# get gene annotation from GENCODE
less gencode.v24.annotation.gtf | grep "\tgene\t" >> ${dir}/gencode_gene_anotation.gtf

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' ${dir}/gencode_gene_anotation.gtf >> ${dir}/gencode_v24_annotation_genes_transcript_id.gtf

# get genomic coordinates of all annotated genes in the human genome in the BED format
gtf2bed < ${dir}/gencode_v24_annotation_genes_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_annotation_genes.bed

# get genomic regions of length ${dist_from_genes} base pairs flanking genes from both sides
# dist_from_genes = 50000
bedtools flank -b ${dist_from_genes} -i ${dir}/gencode_v24_annotation_genes.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_annotation_genes_flanks.bed

# merge regions containing genes and their flanking regions
cat ${dir}/gencode_v24_annotation_genes.bed ${dir}/gencode_v24_annotation_genes_flanks.bed >> ${dir}/gencode_v24_annotation_genes_with_flanks.bed

sortBed -i ${dir}/gencode_v24_annotation_genes_with_flanks.bed >> ${dir}/gencode_v24_annotation_genes_with_flanks_sorted.bed

bedtools merge -i ${dir}/gencode_v24_annotation_genes_with_flanks_sorted.bed >> ${dir}/gencode_v24_annotation_genes_with_flanks_merged.bed

#==================================================================================
#
#				   MicroRNAs
#
#==================================================================================
echo "Distance from microRNAs = ${dist_from_micrornas}" bp
# dist_from_micrornas = 300000
mkdir tmp/micrornas
dir=tmp/micrornas

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
bedtools flank -b 300000 -i ${dir}/ENCFF417SFH_no_chrEBV2.bed -g data/chromInfo_hg38.txt >> ${dir}/ENCFF417SFH_flanks.bed
bedtools flank -b 300000 -i ${dir}/ENCFF850OSS_no_chrEBV2.bed -g data/chromInfo_hg38.txt >> ${dir}/ENCFF850OSS_flanks.bed

# merge regions containing microRNAs and their flanking regions
cat ${dir}/ENCFF850OSS_no_chrEBV.bed ${dir}/ENCFF850OSS_flanks.bed >> ${dir}/ENCFF850OSS_with_flanks.bed
cat ${dir}/ENCFF417SFH_no_chrEBV.bed ${dir}/ENCFF417SFH_flanks.bed >> ${dir}/ENCFF417SFH_with_flanks.bed

sortBed -i ${dir}/ENCFF417SFH_with_flanks.bed >> ${dir}/ENCFF417SFH_with_flanks_sorted.bed
sortBed -i ${dir}/ENCFF850OSS_with_flanks.bed >> ${dir}/ENCFF850OSS_with_flanks_sorted.bed

bedtools merge -i ${dir}/ENCFF417SFH_with_flanks_sorted.bed >> ${dir}/ENCFF417SFH_with_flanks_merged.bed
bedtools merge -i ${dir}/ENCFF850OSS_with_flanks_sorted.bed >> ${dir}/ENCFF850OSS_with_flanks_merged.bed

#==================================================================================
#
#				     tRNAs
#
#==================================================================================
echo "Distance from tRNAs = ${dist_from_trnas}" bp

mkdir tmp/trnas
dir=tmp/trnas

# get tRNA annotation from GENCODE
gtf2bed < gencode.v24.tRNAs.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_tRNAs.bed

# dist_from_trnas = 300000
# get genomic regions of length ${dist_from_trnas} base pairs flanking tRNAs from both sides
bedtools flank -b ${dist_from_trnas} -i ${dir}/gencode_v24_tRNAs.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_tRNAs_flanks.bed

# merge regions containing tRNAs and their flanking regions
cat ${dir}/gencode_v24_tRNAs.bed ${dir}/gencode_v24_tRNAs_flanks.bed >> ${dir}/gencode_v24_tRNAs_with_flanks.bed

sortBed -i ${dir}/gencode_v24_tRNAs_with_flanks.bed >> ${dir}/gencode_v24_tRNAs_with_flanks_sorted.bed

bedtools merge -i ${dir}/gencode_v24_tRNAs_with_flanks_sorted.bed >> ${dir}/gencode_v24_tRNAs_with_flanks_merged.bed

#==================================================================================
#
#			      Long-non-coding RNAs
#
#==================================================================================

echo "Distance from lncRNAs = ${dist_from_lncRNAs}" bp

mkdir tmp/lncrnas
dir=tmp/lncrnas

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v24.long_noncoding_RNAs.gtf >> ${dir}/gencode_v24_long_noncoding_RNAs_transcript_id.gtf

gtf2bed < ${dir}/gencode_v24_long_noncoding_RNAs_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_long_noncoding_RNAs.bed

# dist_from_lncRNAs=300000
bedtools flank -b ${dist_from_lncRNAs} -i ${dir}/gencode_v24_long_noncoding_RNAs.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_long_noncoding_RNAs_flanks.bed

cat ${dir}/gencode_v24_long_noncoding_RNAs.bed ${dir}/gencode_v24_long_noncoding_RNAs_flanks.bed >> ${dir}/gencode.v24.long_noncoding_RNAs_with_flanks.bed

sortBed -i ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks.bed >> ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_sorted.bed

bedtools merge -i ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_sorted.bed >> ${dir}/gencode_v24_long_noncoding_RNAs_with_flanks_merged.bed

#==================================================================================
#
#				   Oncogenes
#
#==================================================================================

echo "Distance from oncogenes = ${dist_from_oncogenes}" bp

mkdir tmp/oncogenes
dir=tmp/oncogenes

# get GENCODE gene annotation for oncogenes from COSMIC (Cancer Gene Census)
grep -f data/Oncogenes_id_list.txt gencode.v24.annotation.gtf >> tmp/oncogenes/gencode_oncogenes_annotation.gtf

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' ${dir}/gencode_oncogenes_annotation.gtf >> ${dir}/gencode_oncogenes_annotation_transcript_id.gtf

# get genomic coordinates of the oncogenes
gtf2bed < ${dir}/gencode_oncogenes_annotation_transcript_id.gtf | awk -v OFS="\t" '{print $1, $2, $3}' >> ${dir}/gencode_v24_annotation_oncogenes.bed

# get genomic regions of length ${dist_from_oncogenes} base pairs flanking oncogenes from both sides
# dist_from_oncogenes = 300000
bedtools flank -b ${dist_from_oncogenes} -i ${dir}/gencode_v24_annotation_oncogenes.bed -g data/chromInfo.txt >> ${dir}/gencode_v24_annotation_oncogenes_flanks.bed

# merge regions containing oncogenes and their flanking regions
cat ${dir}/gencode_v24_annotation_oncogenes.bed ${dir}/gencode_v24_annotation_oncogenes_flanks.bed >> ${dir}/gencode_v24_annotation.oncogenes_with_flanks.bed

sortBed -i ${dir}/gencode_v24_annotation.oncogenes_with_flanks.bed >> ${dir}/gencode_v24_annotation_oncogenes_with_flanks_sorted.bed

bedtools merge -i ${dir}/gencode_v24_annotation_oncogenes_with_flanks_sorted.bed >> ${dir}/gencode_v24_annotation_oncogenes_with_flanks_merged.bed

#==================================================================================
#
#				    Enhancers
#
#==================================================================================

echo "Distance from enhancers = ${dist_from_enhancers}" bp

mkdir tmp/enhancers
dir=tmp/enhancers

# merge regions containing enhancers from the ROADMAP project
for filename in data/Roadmap/*.bed; do
    cat $filename >> ${dir}/All_enhancers.bed
done

sortBed -i ${dir}/All_enhancers.bed >> ${dir}/All_enhancers_sorted.bed

bedtools merge -i ${dir}/All_enhancers_sorted.bed >> ${dir}/All_enhancers_merged.bed

# get genomic regions of length ${dist_from_enhancers} base pairs flanking enhancers from both sides
# dist_from_enhancers=300000
bedtools flank -b ${dist_from_enhancers} -i ${dir}/All_enhancers_merged.bed -g data/chromInfo_hg38.txt >> ${dir}/All_enhancers_flanks.bed

# merge regions containing enhancers and their flanking regions
cat ${dir}/All_enhancers_merged.bed ${dir}/All_enhancers_flanks.bed >> ${dir}/All_enhancers_with_flanks.bed

sortBed -i ${dir}/All_enhancers_with_flanks.bed >> ${dir}/All_enhancers_with_flanks_sorted.bed

bedtools merge -i ${dir}/All_enhancers_with_flanks.bed >> ${dir}/All_enhancers_with_flank_merged.bed

#==================================================================================
#
#		 	Union of all genomic regions to avoid
#
#==================================================================================

