#!/bin/bash

GENE_SYMBOL_FILE="/Users/hansep/PycharmProjects/diachrscripts/rscripts/JAV_MK_RALT_0.0019_interactions_with_genesymbols.dc.tsv"
OUT_PREFIX="MK_RALT_0.0019"

#GENE_SYMBOL_FILE="/Users/hansep/PycharmProjects/diachrscripts/rscripts/JAV_ERY_RALT_0.0018_interactions_with_genesymbols.dc.tsv"
#OUT_PREFIX="ERY_RALT_0.0018"
#/home/peter/storage_1/VPV_data/hg38/hg38.fa

# MEME - ERY_RALT_0.0018 - STAA vs. URAA - same strand - repeats masked
# MEME - ERY_RALT_0.0018 - STAA vs. URAA - same strand

# DREME - ERY_RALT_0.0018 - STAA vs. URAA - same strand - repeats masked
# DREME - ERY_RALT_0.0018 - STAA vs. URAA - same strand

BED_TOOLS="bedtools"
GENOME_SEQ="/Users/hansep/data/hg38/hg38.fa"

echo "Extracting directed regions ..."
awk '{if(($3=="S" || $3=="T") && $6=="AA"  && ($8=="+/+" || $8=="-/-")){split($1,A,";");split(A[1],D1,":");split(A[2],D2,":");split(D1[2],D1C,"-");split(D2[2],D2C,"-");print D1[1]"\t"D1C[1]"\t"D1C[2]"\tD1|"$3"|"$2"|"$8"|"$7"\n"D2[1]"\t"D2C[1]"\t"D2C[2]"\tD2|"$3"|"$2"|"$8"|"$7}}' $GENE_SYMBOL_FILE > $OUT_PREFIX\_DAA.bed
echo "... done."

echo "Extracting undirected regions ..."
awk '{if($3=="URAA" && $6=="AA"  && ($8=="+/+" || $8=="-/-")){split($1,A,";");split(A[1],D1,":");split(A[2],D2,":");split(D1[2],D1C,"-");split(D2[2],D2C,"-");print D1[1]"\t"D1C[1]"\t"D1C[2]"\tD1|"$3"|"$2"|"$8"|"$7"\n"D2[1]"\t"D2C[1]"\t"D2C[2]"\tD2|"$3"|"$2"|"$8"|"$7}}' $GENE_SYMBOL_FILE > $OUT_PREFIX\_URAA.bed
echo "... done."

echo "Coverting BED to FASTA format for directed interactions ..."
$BED_TOOLS getfasta -name -fi $GENOME_SEQ -bed $OUT_PREFIX\_DAA.bed > $OUT_PREFIX\_DAA.fasta
echo "... done."

echo "Coverting BED to FASTA format for undirected interactions ..."
$BED_TOOLS getfasta -name -fi $GENOME_SEQ -bed $OUT_PREFIX\_URAA.bed > $OUT_PREFIX\_URAA.fasta
echo "... done."

echo "Masking repetitive sequences for directed interactions ..."
awk '{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}' $OUT_PREFIX\_DAA.fasta >  $OUT_PREFIX\_DAA_masked.fasta
echo "... done."

echo "Masking repetitive sequences for undirected interactions ..."
awk '{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}' $OUT_PREFIX\_URAA.fasta >  $OUT_PREFIX\_URAA_masked.fasta
echo "... done."