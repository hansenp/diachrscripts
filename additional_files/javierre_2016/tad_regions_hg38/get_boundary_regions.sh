#!/bin/bash

EXT_SIZE=10000

MERGE_THRESHOLD=5000

tad_region_files_list=(\
"hglft_genome_TADs_MK_hg38.bed" \
"hglft_genome_TADs_ERY_hg38.bed" \
"hglft_genome_TADs_NEU_hg38.bed" \
"hglft_genome_TADs_MON_hg38.bed" \
"hglft_genome_TADs_MAC_M0_hg38.bed" \
"hglft_genome_TADs_NB_hg38.bed" \
"hglft_genome_TADs_NCD4_hg38.bed" \
"hglft_genome_TADs_NCD8_hg38.bed")


echo "track name=\"All TAD boundaries\" description=\"All TAD boundaries\" visibility=4" > all_tad_boundaries.bed
for file in "${tad_region_files_list[@]}"
do
  tail -n +2 $file | awk '{print $1"\t"$2-'${EXT_SIZE}'"\t"$2+'${EXT_SIZE}'"\n"$1"\t"$3-'${EXT_SIZE}'"\t"$3+'${EXT_SIZE}'}' >> all_tad_boundaries.bed
done


if [ ! -e 'sorted_tad_boundaries.bed' ]; then
  sort -k1,1 -k2,2n all_tad_boundaries.bed > sorted_tad_boundaries.bed
fi

bedtools merge -i sorted_tad_boundaries.bed -d ${MERGE_THRESHOLD} > merged_tad_boundaries.bed