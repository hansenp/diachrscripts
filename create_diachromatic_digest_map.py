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

print("\t[INFO] Read " + str(len(enriched_digests_set)) + " digests ...")
print("[INFO] ... done.")


### Iterate input DigestMap and write new DigestMap with enrichment status flags overwritten
############################################################################################

#with open(diachromatic_digest_map, 'rt') as fp:
    #for line in fp:
        #pass
