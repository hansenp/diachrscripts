#!/usr/bin/env python

"""
This script takes a BED file with interaction coordinates that consist of the start coordinate of the first and the end
coordinate of the second interacting digest. The coordinates are added to a set and enhanced interaction file is then
iterated, whereby all interactions that are not contained in the set will be discarded.
"""


import argparse
import gzip
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Takes a BED file with interaction coordinates and filters an enhanced interaction file for these interactions.')

parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--bed-interaction-file', required=True, help='BED file with start coordinates of the first and end coordinate of the second interacting digest.')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
bed_interaction_file = args.bed_interaction_file
enhanced_interaction_file = args.enhanced_interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --bed-interaction-file: " + bed_interaction_file)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)

enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_filtered.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')

### Start execution
###################

print("[INFO] Iterating BED file ...")
n_interaction_chicago = 0
valid_interaction_set = set()
with open(bed_interaction_file, 'r' + 't') as fp:

    for line in fp:

        # Count total number of interactions
        n_interaction_chicago += 1

        # Report progress
        if n_interaction_chicago % 10000 == 0:
            print("\t\t[INFO]", n_interaction_chicago, "interactions processed ...")

        # Create key from coordinates and add key to set
        chr, sta, end = line.split('\t')
        key = (chr + ":" + sta + "-" + end).rstrip()
        valid_interaction_set.add(key)

fp.close()

print("[INFO] Iterating enhanced interaction file ...")
cnt_selected = 0
n_interaction_diachromatic = 0
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    for line in fp:

        # Count total number of Diachromatic interactions
        n_interaction_diachromatic += 1

        # Report progress
        if n_interaction_diachromatic % 100000 == 0:
            print("\t\t[INFO]", n_interaction_diachromatic, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Create key from coordinates and check whether the key is contained in the set
        key = chr_a + ":" + str(sta_a) + "-" + str(end_b)

        if key in valid_interaction_set:
            enhanced_interaction_stream_output.write(line.rstrip() + '\n')
            cnt_selected += 1

fp.close()

enhanced_interaction_stream_output.close()

print("Found " + str(cnt_selected) + " out of " + str(n_interaction_chicago) + " CHiCAGO interactions (" + str(cnt_selected/n_interaction_chicago) + ") amongst " + str(n_interaction_diachromatic) + " Diachromatic interactions.")
print("Selected " + str(cnt_selected) + " out of " + str(n_interaction_diachromatic) + " Diachromatic interactions (" + str(cnt_selected/n_interaction_diachromatic) + ").")

print(str(n_interaction_chicago) + "\t" + str(n_interaction_diachromatic) + "\t" + str(cnt_selected))

print("[INFO] ... done.")
