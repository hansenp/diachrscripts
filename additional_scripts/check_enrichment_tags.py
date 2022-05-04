#!/usr/bin/env python

"""
This script verifies that the pass-through of enrichment tags ('E' and 'N') works as intended.
It takes a BED file with coordinates of digests selected for enrichment as well as a Diachromatic interaction
file. A set 'A' is created from the enriched digests in the BED file.
Then the Diachromatic interaction file is iterated, adding digests with an 'E' to a second set ('B') and digests with
an 'N' to a third set 'C'.
If everything works, then the intersection between set 'A' and 'B' must become large when a lot of interactions are read.
Furthermore, 'B' must be completely contained in 'A', and the intersection between A and C must remain empty,
otherwise something is wrong.
"""

import argparse
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


# Parse command line
####################

parser = argparse.ArgumentParser(description='Overwrite enrichment status in column 12 of a Diachromatic digest file.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.', required=True)
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
enriched_digests_file = args.enriched_digests_file
diachromatic_interaction_file = args.diachromatic_interaction_file


print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enriched-digests-file: " + str(enriched_digests_file))
print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)


# Read BED file with coordinates of enriched digests
####################################################

print("[INFO] Reading list with digests selected for enrichment ...")

enriched_digests_set = set()
with open(enriched_digests_file, 'rt') as fp:
    skipped_first_line = fp.readline()
    for line in fp:
        chrom, sta, end = line.rstrip().split('\t')
        enriched_digests_set.add(chrom + '\t' + str(sta) + '\t' + str(end))

n_enriched_digests = len(enriched_digests_set)
print("\t[INFO] Read " + str(n_enriched_digests) + " digests ...")
print("[INFO] ... done.")


# Iterate Diachromatic interaction file and create sets of enriched and not enriched digests
############################################################################################

# Init sets for enriched and not enriched interactions
enriched_digests_interactions = set()
not_enriched_digests_interactions = set()

# Read Diachromatic interaction file
print("[INFO] Reading Diachromatic interaction file")
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file)
print("[INFO] ... done.")

# Iterate interaction objects and create sets of enriched and not enriched digests
print("[INFO] Iterating interaction objects")
n_interaction = 0
for i in interaction_set.interaction_list:

    n_interaction += 1

    if n_interaction % 1000000 == 0:
        print("\t[INFO] Processed " + str(n_interaction) + " interactions ...")

    enrichment_tag_A = i.enrichment_status_tag_pair[0]
    key_A = i.chrA + "\t" + str(i.fromA) + "\t" + str(i.toA)
    if enrichment_tag_A == 'E':
        enriched_digests_interactions.add(key_A)
    else:
        not_enriched_digests_interactions.add(key_A)

    enrichment_tag_B = i.enrichment_status_tag_pair[1]
    key_B = i.chrB + "\t" + str(i.fromB) + "\t" + str(i.toB)
    if enrichment_tag_B == 'E':
        enriched_digests_interactions.add(key_B)
    else:
        not_enriched_digests_interactions.add(key_B)

print("[INFO] ... done.")

print("\nNumber of enriched digests in original set A: " + str(len(enriched_digests_set)))
print("Number of enriched digests in set B derived from interactions ('E'): " + str(len(enriched_digests_interactions)))
print("Number of not enriched digests in set C derived from interactions ('N'): " + str(len(not_enriched_digests_interactions)))

a_intersect_b = enriched_digests_set & enriched_digests_interactions
print("\nIntersect between set A and B (Should be large.): " + str(len(a_intersect_b)))

set_diff_b_minus_a = enriched_digests_interactions -  enriched_digests_set
print("\nSet difference B minus A (Must be empty!): " + str(len(set_diff_b_minus_a)))
if(0 < len(set_diff_b_minus_a)):
    print("Something is wrong! B must be completely contained in A.")

a_intersect_c = enriched_digests_set & not_enriched_digests_interactions
print("\nIntersect between set A and C (Must be empty!): " + str(len(a_intersect_c)))
if(0 < len(a_intersect_c)):
    print("Something is wrong! A and C must not overlap.")
