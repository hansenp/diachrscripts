#!/usr/bin/env python

"""
This script takes a BED file with coordinates of digests selected for enrichment as well as a Diachromatic DigestMap
and creates a new DigestsMap in which digest regions that are contained in the BED file are marked with a '1'
and all other digests are marked with '0'.
"""

import argparse


### Parse command line
######################

parser = argparse.ArgumentParser(description='Overwrite enrichment status in column 12 of a Diachromatic DigestMap.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.', required=True)
parser.add_argument('-i','--diachromatic-digest-map', help='Diachromatic DigestMap created with GOPHER.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
enriched_digests_file = args.enriched_digests_file
diachromatic_digest_map = args.diachromatic_digest_map

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enriched-digests-file: " + str(enriched_digests_file))
print("\t[INFO] --diachromatic-digest-map: " + diachromatic_digest_map)


### Read BED file with coordinates of enriched digests
######################################################

print("[INFO] Reading list with digests selected for enrichment ...")

enriched_digests_set = set()
with open(enriched_digests_file, 'rt') as fp:
    for line in fp:
        chr, sta, end = line.rstrip().split('\t')
        enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))

n_enriched_digests = len(enriched_digests_set)
print("\t[INFO] Read " + str(n_enriched_digests) + " digests ...")
print("[INFO] ... done.")


### Iterate input DigestMap and write new DigestMap with enrichment status flags overwritten
############################################################################################

print("[INFO] Iterating input digest map and writing new digest map ...")

# Create file that contains all the information
out_diachromatic_digest_map_fn =  out_prefix + "_diachromatic_digest_map.txt"
out_diachromatic_digest_map_fh = open(out_diachromatic_digest_map_fn, 'wt')

enriched_digests_found_set = set()
with open(diachromatic_digest_map, 'rt') as fp:
    for line in fp:

        # Skip first line

        # Parse one line of a Diachromatic digest map
        fields = line.rstrip('\n').split('\t')
        chr = fields[0]        # Chromosome
        sta = fields[1]        # Fragment_Start_Position
        end = fields[2]        # Fragment_End_Position
        frN = fields[3]         # Fragment_Number
        rS5 = fields[4]         # 5'_Restriction_Site
        rS3 = fields[5]         # 3'_Restriction_Site
        Len = fields[6]         # Length
        gc5 = fields[7]         # 5'_GC_Content
        gc3 = fields[8]         # 3'_GC_Content
        rp5 = fields[9]         # 5'_Repeat_Content
        rp3 = fields[10]        # 3'_Repeat_Content
        selected = fields[11]  # Selected
        pr5 = fields[12]       # 5'_Probes
        pr3 = fields[13]       # 3'_Probes

        # Check whether digest is in enriched digest set
        SELECTED = 'F'
        PR5 = 0
        PR3 = 0
        key = chr + '\t' + sta + '\t' + end
        if key in enriched_digests_set:
            enriched_digests_found_set.add(key)
            SELECTED = 'T'
            PR5 = 1
            PR3 = 1

        # Write modified line to new digest map
        out_diachromatic_digest_map_fh.write(
            chr + '\t' +
            sta + '\t' +
            end + '\t' +
            frN + '\t' +
            rS5 + '\t' +
            rS3 + '\t' +
            Len + '\t' +
            gc5 + '\t' +
            gc3 + '\t' +
            rp5 + '\t' +
            rp3 + '\t' +
            SELECTED + '\t' +
            str(PR5) + '\t' +
            str(PR3)
            )

out_diachromatic_digest_map_fh.close()

n_enriched_digests_matched = len(enriched_digests_found_set)
n_enriched_digests_unmatched = n_enriched_digests - n_enriched_digests_matched

print("\t[INFO] From the " + str(n_enriched_digests) + " enriched digests in the input file, "
      + str(n_enriched_digests_matched) + " were found in the digest map.")
print("\t[INFO] " + str(n_enriched_digests_unmatched) + " digest were not found.")

print("[INFO] ... done.")

print("[INFO] Trying to get more information about digests that were not found in the digest map ...")

add_range = 100
id_not_found = 0
for d in enriched_digests_set:
    if d not in enriched_digests_found_set:
        id_not_found += 1
        print("------------------------------------------------------")
        print("DIGEST NOT FOUND " + str(id_not_found) + '\n')

        print("Digest coordinates:")
        print(d + '\n')

        CHR = d.split('\t')[0]
        STA = int(d.split('\t')[1])
        END = int(d.split('\t')[2])
        LEN = END - STA

        print("Digest length:")
        print(str(LEN) + '\n')

        with open(diachromatic_digest_map, 'rt') as fp:
            for line in fp:

                fields = line.rstrip('\n').split('\t')
                chr = fields[0]        # Chromosome
                sta = fields[1]        # Fragment_Start_Position
                end = fields[2]        # Fragment_End_Position

                if chr == CHR and STA - add_range <= int(sta) and int(end) <= END + add_range:
                    print(line)

print("[INFO] ... done.")
# chr15	23033835	23052286
# chr15	23033832	23052283	2142	HindIII	HindIII	18452	0.005	0.005	0.013	0.014	F	0	0


# --out-prefix
# results/tmp_results/create_digest_map/JAV_HindIII_hg19
# --enriched-digests-file
# data/tmp_example_input_data/digest_map/hg19/Digest_Human_HindIII_baits_e75_ID.baitmap.bed
# --diachromatic-digest-map
# data/tmp_example_input_data/digest_map/hg19/digest_map_hg19_hg19_DigestedGenome.txt


# --out-prefix
# results/tmp_results/create_digest_map/JAV_HindIII_hg19
# --enriched-digests-file
# additional_files/javierre_2016/baited_digest_regions/Digest_Human_HindIII_baits_e75_ID.baitmap.hg38.bed
# --diachromatic-digest-map
# data/tmp_example_input_data/digest_map/hg38/digest_map_hg38_hg38_DigestedGenome.txt