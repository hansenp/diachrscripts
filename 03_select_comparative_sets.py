#!/usr/bin/env python

"""
This script takes an enhanced interaction file supplemented with digest associated gene symbols and TSS as well as
directionality P-values. The script iterates the interaction file three times:

   1. In a first pass, the set of directed interactions (DI) is defined based on a specified P-value threshold.

   2. In a second pass, bad digests that interact with digests involved in DI (DD) are identified.

   3. In a third pass, the set of exclusive undirected interactions (EUI) and associated digests (EUD) is identified
      and DI and EUI are written to a enhanced interaction formatted file.

Finally, summary statistics about interaction and digest sets are reported.
"""


import argparse
import math
import gzip
import diachrscripts_toolkit


### Parse command line
######################

parser = argparse.ArgumentParser(description='Select comparative sets of interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--p-value-threshold', help='P-value threshold for directed interactions.', default=0.01)
parser.add_argument('--min-digest-dist', help='All interactions with smaller distances will be discarded.', default=20000)
args = parser.parse_args()

out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
p_value_threshold = float(args.p_value_threshold)
min_digest_dist = int(args.min_digest_dist)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --p-value-cutoff: " + str(p_value_threshold))
print("\t[INFO] --min-digest-dist: " + str(min_digest_dist))


### Define auxiliary functions
##############################


### Prepare variables, data structures and streams for output files
###################################################################

# P-value threshold
neg_log_p_val_thresh =  round((-1) * math.log(p_value_threshold), 2)

# number of directed interactions and involved digest
dir_inter_num = 0
dir_inter_involved_dig_num = None

# number of undirected interactions and involved digest
undir_inter_num = 0
undir_inter_involved_dig_num = None

# number of indefinable interactions
indef_inter_num = 0

# number of exclusive undirected interactions
undir_exc_inter_num = 0

# Set containing all digests involved in directed interactions
dir_dig_set = set()

# Set containing all digests involved in undirected interactions
undir_dig_set = set()

# Set containing all digests involved in exclusive undirected interactions
undir_exc_dig_set = set()

# Set containing digests of undirected interactions that interact with digests involved in directed interactions
bad_dig_set = set()

# Smallest n that required for significance given the P-value threshold
smallest_n = diachrscripts_toolkit.find_indefinable_n(p_value_threshold, verbose = False)


### 1st pass: Identifiy directed interactions (DI) based on P-value threshold
#############################################################################

print("[INFO] 1st pass: Identify directed interactions (DI) based on P-value threshold ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Directed interactions
        if neg_log_p_val_thresh <= neg_log_p_value:
            dir_inter_num += 1
            dir_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            dir_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")


### 2nd pass: Identifiy bad digests that interact with digests involved in DI
#############################################################################

print("[INFO] 2nd pass: Identifiy bad digests that interact with digests involved in DI ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Skip and count indefinable interactions
        if rp_total < smallest_n:
            indef_inter_num +=1
            line = fp.readline()
            continue

        # Skip directed interactions. All remaining interactions must be undirected.
        if neg_log_p_val_thresh <= neg_log_p_value:
            line = fp.readline()
            continue

        undir_inter_num += 1

        d1_coords = chr_a + "\t" + str(sta_a) + "\t" + str(end_a)
        d2_coords = chr_b + "\t" + str(sta_b) + "\t" + str(end_b)

        # Add digests involved in undirected interactions to set
        undir_dig_set.add(d1_coords)
        undir_dig_set.add(d2_coords)

        # Check if the first digest of the interaction interacts with a digest involved in DI
        if d2_coords in dir_dig_set:
            bad_dig_set.add(d1_coords)

        # Check if the second digest of the interaction interacts with a digest involved in DI
        if d1_coords in dir_dig_set:
            bad_dig_set.add(d2_coords)

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")


### 3rd: Identifiy exclusive undirected interctions and associated digests
##########################################################################

print("[INFO] 3rd: Identifiy exclusive undirected interctions and associated digests ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Skip and count indefinable interactions
        if rp_total < smallest_n:
            indef_inter_num +=1
            line = fp.readline()
            continue

        # Skip directed interactions. All remaining interactions must be undirected.
        if neg_log_p_val_thresh <= neg_log_p_value:
            line = fp.readline()
            # print line with
            continue

        d1_coords = chr_a + "\t" + str(sta_a) + "\t" + str(end_a)
        d2_coords = chr_b + "\t" + str(sta_b) + "\t" + str(end_b)

        #if not(d1_coords in dir_dig_set or d1_coords in bad_dig_set or d2_coords in dir_dig_set or d2_coords in bad_dig_set):
        if not (d1_coords in dir_dig_set or d2_coords in dir_dig_set):
            undir_exc_dig_set.add(d1_coords)
            undir_exc_dig_set.add(d2_coords)
            undir_exc_inter_num += 1

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")


### Finish: Output statistics about interaction and digest sets
###############################################################

print("[INFO] Finish: Output statistics about interaction and digest sets ...")

print()

print("\t[INFO] Number of directed interactions: " + str(dir_inter_num))
dir_inter_involved_dig_num = len(dir_dig_set)
print("\t[INFO] Number of digests that are involved in directed interactions: " + str(dir_inter_involved_dig_num))
dir_inter_connectivity = "{0:.2f}".format(1 - (dir_inter_involved_dig_num / (2 * dir_inter_num)))
print("\t[INFO] Connectivity factor for directed interactions: " + str(dir_inter_connectivity))

print()

print("\t[INFO] Number of indefinable interactions: " + str(indef_inter_num))
print("\t[INFO] Number of undirected interactions: " + str(undir_inter_num))
undir_inter_involved_dig_num = len(undir_dig_set)
print("\t[INFO] Number of digests that are involved in undirected interactions: " + str(undir_inter_involved_dig_num))
undir_inter_connectivity = "{0:.2f}".format(1 - (undir_inter_involved_dig_num / (2 * undir_inter_num)))
print("\t[INFO] Connectivity factor for undirected interactions: " + str(undir_inter_connectivity))

print()

print("\t[INFO] Number of bad digests: " + str(len(bad_dig_set)))

print()

print("\t[INFO] Number of undirected interactions: " + str(undir_inter_num))
print("\t[INFO] Number of exclusive undirected interactions: " + str(undir_exc_inter_num))
undir_exc_inter_involved_dig_num = len(undir_exc_dig_set)
print("\t[INFO] Number of digests involved exclusive undirected interactions: " + str(undir_exc_inter_involved_dig_num))
undir_exc_inter_connectivity = "{0:.2f}".format(1 - (undir_exc_inter_involved_dig_num / (2 * undir_exc_inter_num)))
print("\t[INFO] Connectivity factor for exclusive undirected interactions: " + str(undir_exc_inter_connectivity))







