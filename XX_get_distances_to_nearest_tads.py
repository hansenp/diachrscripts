#!/usr/bin/env python

"""
This script takes as input an enhanced interaction file, a BED file with TAD regions and, optionally, a file with the
coordinates of digests that will be excluded from the analysis. Three digest sets are derived
from the interaction file:

   1. Digests that are involved in directed interactions only (DD).
   2. Digests that are involved in undirected interactions only (AD).
   3. Digests that are involved in directed and undirected interactions (UD).

The TAD regions are stored in an object of the class 'TADs' that is essentially a hash map. The names of the
chromosomes serve as keys and for each key there is a sorted list of boundary coordinates. The class has a function that
returns the coordinate of the next TAD boundary for a given coordinate. This is accomplished through binary search.

For each digest, the distance to the nearest TAD boundary is determined and added to a list that corresponds to one of
the three digest sets (DD, AD and UD). The list with distances is written to a TSV file that can be used for further
analyses in R (histograms). The digest regions for the three set are written to separate BED files and can be used for
downnstream analyses. If a digest is on a chromosome that does not have a TAD region, no distance is taken into account
(TSV file), but the digest is still written to the BED file. For the capture Hi-C dataset of Javierre et al. 2016,
this only occurs for the MK dataset. 47 digest on chromosome Y are affected.

Black list: We used the '.baitmap' file for CHiCAGO used in the original publication (Javierre 2016) to exclude digests
that were selected for enrichment. We downloaded the file from here:

https://osf.io/u8tzp/wiki/home/

We used UCSC's liftOver tool to convert the coordinates of baited digests from hg19 to hg38.

"""


import argparse
import gzip
import diachrscripts_toolkit as dclass
from bisect import bisect_left
import numpy as np


### Parse command line
######################

parser = argparse.ArgumentParser(description='Get distances to nearest TAD boundary for diegstes involved in directed interactions and undirected reference interactions not involved in directed interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--tad-region-bed-file', help='BED file with TAD regions.', required=True)
parser.add_argument('--digest-black-list', help='BED file with digests to be excluded.')

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
tad_region_bed_file = args.tad_region_bed_file
digest_black_list = args.digest_black_list

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --tad-region-bed-file: " + tad_region_bed_file)
if digest_black_list != None:
    print("\t[INFO] --digest-black-list: " + digest_black_list)


### Definition of the class 'Tad_boundaries'
############################################

class Tad_boundaries:

    # Class to store TAD boundaries with a function that returns the next boundary for a given coordinate

    # Attributes

    n_tads = 0   # Total number of TAD boundaries

    chr_dict = {}   # Dictionary that has a list of boundary coordinates for each chromosome

    # Initializer

    def __init__(self, tad_region_bed_file):

        with open(tad_region_bed_file, 'rt') as fp:
            line = fp.readline().rstrip()
            while line:

                # Parse line
                chr_key, boundary_1, boundary_2 = line.split('\t')

                # Count TAD boundaries
                self.n_tads = self.n_tads + 2

                # Add boundary coordinates to dictionary
                if chr_key in self.chr_dict:
                    self.chr_dict[chr_key].append(int(boundary_1))
                    self.chr_dict[chr_key].append(int(boundary_2))
                else:
                    self.chr_dict[chr_key] = [int(boundary_1), int(boundary_2)]

                # Go to next line
                line = fp.readline().rstrip()

        fp.close()

        # Sort lists with coordinates for each chromosome
        for chr_key in self.chr_dict:
            self.chr_dict[chr_key] = sorted(self.chr_dict[chr_key])

    def get_nearest_tad_boundary(self, chr_key, coord):

        # Get the index of x (if in the list) or the index of next larger coordinate
        boundary_to_the_right_idx = bisect_left(self.chr_dict[chr_key],int(coord))

        if len(self.chr_dict[chr_key]) <= boundary_to_the_right_idx:
            # There is no boundary to the right
            return self.chr_dict[chr_key][boundary_to_the_right_idx - 1]

        if boundary_to_the_right_idx == 0:
            # There is no boundary to the left
            return self.chr_dict[chr_key][boundary_to_the_right_idx]

        if self.chr_dict[chr_key][boundary_to_the_right_idx] == coord:
            # coord is a TAD boundary and the distance is zero
            return coord
        else:
            # coord is between two TAD boundaries, one to the left and one to the right
            boundary_to_the_left_idx = boundary_to_the_right_idx - 1

            if (self.chr_dict[chr_key][boundary_to_the_right_idx] - coord) < (coord - self.chr_dict[chr_key][boundary_to_the_left_idx]):
                # The boundary to the right is closer
                return self.chr_dict[chr_key][boundary_to_the_right_idx]
            else:
                # The boundary to the left is closer
                return self.chr_dict[chr_key][boundary_to_the_left_idx]

    def get_distance_to_nearest_tad_boundary(self, chr_key, coord):

        # Check whether there are TAD boundaries on this chromosome
        if chr_key not in self.chr_dict:
            return -1

        # Determine distance to next TAD boundary
        nearest_tad_boundary_coord = self.get_nearest_tad_boundary(chr_key,coord)
        if coord <= nearest_tad_boundary_coord:
            return nearest_tad_boundary_coord - coord
        else:
            return coord - nearest_tad_boundary_coord


### Begin execution
###################

# Init object for TAD boundaries
TAD_boundaries = Tad_boundaries(tad_region_bed_file)

# Init digest sets
digests_from_dir_inter = set()
digests_from_undir_inter = set()
digest_to_be_excluded = set()

# Read digest regions to be excluded from file
if digest_black_list != None:
    with open(digest_black_list, 'rt') as fp:
        line = fp.readline().rstrip()
        while line:

            # Parse line
            chr, sta, end = line.split('\t')

            # Add to set
            digest_to_be_excluded.add(chr + '\t' + str(sta) + '\t' + str(end))

            # Go to next line
            line = fp.readline().rstrip()


# Determine comparative digest sets
print("[INFO] Determining comparative digest sets ...")
print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Add digests to set for directed interactions
        if interaction_category == "DIII" or interaction_category == "DIAI" or interaction_category == "DIAA":

            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) not in digest_to_be_excluded:
                digests_from_dir_inter.add(chr_a + '\t' + str(sta_a) + '\t' + str(end_a))
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) not in digest_to_be_excluded:
                digests_from_dir_inter.add(chr_b + '\t' + str(sta_b) + '\t' + str(end_b))

        # Add digests to set for undirected interactions
        if interaction_category == "UIRII" or interaction_category == "UIRAI" or interaction_category == "UIRAA":
            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) not in digest_to_be_excluded:
                digests_from_undir_inter.add(chr_a + '\t' + str(sta_a) + '\t' + str(end_a))
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) not in digest_to_be_excluded:
                digests_from_undir_inter.add(chr_b + '\t' + str(sta_b) + '\t' + str(end_b))

        line = fp.readline()

fp.close()
print("\t[INFO] ... done.")

dd_digest_set = digests_from_dir_inter - digests_from_undir_inter
ud_digest_set = digests_from_undir_inter - digests_from_dir_inter
ad_digest_set = digests_from_dir_inter & digests_from_undir_inter

print("\t[INFO] Number of digests involved in directed interactions: " + str(len(digests_from_dir_inter)))
print("\t[INFO] Number of digests involved in undirected interactions: " + str(len(digests_from_undir_inter)))
print("\t[INFO] Number of digests involved in directed interactions only (DD): " + str(len(dd_digest_set)))
print("\t[INFO] Number of digests involved in directed and undirected interactions (AD): " + str(len(ad_digest_set)))
print("\t[INFO] Number of digests involved in undirected interactions only (UD): " + str(len(ud_digest_set)))

print("[INFO] ... done.")

print("[INFO] Determining UNIQUE distances between digests and next TAD ...")

print("\t[INFO] For digests involved in directed interactions only ...")
total_dist_dd = 0
dd_dist_array = []
for digest in dd_digest_set:
    chr_key, sta, end = digest.split('\t')
    digest_center_pos = int(sta) + int((int(end) - int(sta))/2)
    dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_key,digest_center_pos)
    if dist == -1:
        print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
    else:
        total_dist_dd = total_dist_dd + dist
        dd_dist_array.append(dist)

print("\t[INFO] For digests involved in undirected interactions only ...")
total_dist_ud= 0
ud_dist_array = []
for digest in ud_digest_set:
    chr_key, sta, end = digest.split('\t')
    digest_center_pos = int(sta) + int((int(end) - int(sta))/2)
    dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_key,digest_center_pos)
    if dist == -1:
        print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
    else:
        total_dist_ud = total_dist_ud + dist
        ud_dist_array.append(dist)

print("\t[INFO] For digests involved in undirected interactions only ...")
total_dist_ad = 0
ad_dist_array = []
for digest in ad_digest_set:
    chr_key, sta, end = digest.split('\t')
    digest_center_pos = int(sta) + int((int(end) - int(sta))/2)
    dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_key,digest_center_pos)
    if dist == -1:
        print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
    else:
        total_dist_ad = total_dist_ad + dist
        ad_dist_array.append(dist)

print("\t[INFO] Median distance between digests involved in directed interactions only and next TAD boundary: " + str(np.median(dd_dist_array)))
print("\t[INFO] Median distance between digests involved in undirected interactions only and next TAD boundary: " + str(np.median(ud_dist_array)))
print("\t[INFO] Median distance between digests involved in directed and undirected interactions and next TAD boundary: " + str(np.median(ad_dist_array)))

print("\t[INFO] Writing distances to files ...")

tab_stream_dist_to_tad_dd_output = open(out_prefix + "_dist_to_tad_dd_unique.tsv", 'wt')
for dist in dd_dist_array:
    tab_stream_dist_to_tad_dd_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_dd_output.close()

tab_stream_dist_to_tad_ud_output = open(out_prefix + "_dist_to_tad_ud_unique.tsv", 'wt')
for dist in ud_dist_array:
    tab_stream_dist_to_tad_ud_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_ud_output.close()

tab_stream_dist_to_tad_ad_output = open(out_prefix + "_dist_to_tad_ad_unique.tsv", 'wt')
for dist in ad_dist_array:
    tab_stream_dist_to_tad_ad_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_ad_output.close()

print("[INFO] ... done.")

print("[INFO] Writing BED files containing digest regions involved in directed or undirected interactions only ...")

bed_stream_dist_to_tad_dd_output = open(out_prefix + "_dist_to_tad_dd_unique.bed", 'wt')
for digest_region in dd_digest_set:
    bed_stream_dist_to_tad_dd_output.write(digest_region + '\n')
bed_stream_dist_to_tad_dd_output.close()

bed_stream_dist_to_tad_ud_output = open(out_prefix + "_dist_to_tad_ud_unique.bed", 'wt')
for digest_region in ud_digest_set:
    bed_stream_dist_to_tad_ud_output.write(digest_region + '\n')
bed_stream_dist_to_tad_ud_output.close()

bed_stream_dist_to_tad_ad_output = open(out_prefix + "_dist_to_tad_ad_unique.bed", 'wt')
for digest_region in ad_digest_set:
    bed_stream_dist_to_tad_ad_output.write(digest_region + '\n')
bed_stream_dist_to_tad_ad_output.close()

print("[INFO] ... done.")


print("[INFO] Determining REDUNDANT distances between digests and next TAD ...")
dd_dist_array_redundant = []
ud_dist_array_redundant = []
ad_dist_array_redundant = []
dd_interaction_end_array_redundant = []
ud_interaction_end_array_redundant = []
ad_interaction_end_array_redundant = []
n_dd_interaction_ends = 0
n_ud_interaction_ends = 0
n_ad_interaction_ends = 0
print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        digest_coord_key_a = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        digest_coord_key_b = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        # Add distances to list of distances between directed interaction ends and next TAD boundaries
        if interaction_category == "DIII" or interaction_category == "DIAI"  or interaction_category == "DIAA":

            # First digest of interaction
            if digest_coord_key_a in dd_digest_set:
                digest_center_pos = int(sta_a) + int((int(end_a) - int(sta_a)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_a, digest_center_pos)
                if dist != -1:
                    dd_dist_array_redundant.append(dist)
                    dd_interaction_end_array_redundant.append(digest_coord_key_a)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_a + "!")
                n_dd_interaction_ends += 1

            # Second digest of interaction
            if digest_coord_key_b in dd_digest_set:
                digest_center_pos = int(sta_b) + int((int(end_b) - int(sta_b)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_b, digest_center_pos)
                if dist != -1:
                    dd_dist_array_redundant.append(dist)
                    dd_interaction_end_array_redundant.append(digest_coord_key_b)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_b + "!")
                n_dd_interaction_ends += 1


        # Add distances to list of distances between undirected interaction ends and next TAD boundaries
        if interaction_category == "UIRII" or interaction_category == "UIRAI" or interaction_category == "UIRAA":

            # First digest of interaction
            if digest_coord_key_a in ud_digest_set:
                digest_center_pos = int(sta_a) + int((int(end_a) - int(sta_a)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_a, digest_center_pos)
                if dist != -1:
                    ud_dist_array_redundant.append(dist)
                    ud_interaction_end_array_redundant.append(digest_coord_key_a)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
                n_ud_interaction_ends += 1

            # Second digest of interaction
            if digest_coord_key_b in ud_digest_set:
                digest_center_pos = int(sta_b) + int((int(end_b) - int(sta_b)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_b, digest_center_pos)
                if dist != -1:
                    ud_dist_array_redundant.append(dist)
                    ud_interaction_end_array_redundant.append(digest_coord_key_b)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
                n_ud_interaction_ends += 1

        # Add distances to list of distances between ambiguous interaction ends and next TAD boundaries
        if interaction_category == "DIII" or interaction_category == "DIAI"  or interaction_category == "DIAA":

            # First digest of interaction
            if digest_coord_key_a in ad_digest_set:
                digest_center_pos = int(sta_a) + int((int(end_a) - int(sta_a)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_a, digest_center_pos)
                if dist != -1:
                    ad_dist_array_redundant.append(dist)
                    ad_interaction_end_array_redundant.append(digest_coord_key_a)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
                n_ad_interaction_ends += 1

            # Second digest of interaction
            if digest_coord_key_b in ad_digest_set:
                digest_center_pos = int(sta_b) + int((int(end_b) - int(sta_b)) / 2)
                dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_b, digest_center_pos)
                if dist != -1:
                    ad_dist_array_redundant.append(dist)
                    ad_interaction_end_array_redundant.append(digest_coord_key_b)
                else:
                    print("\t\t[Warning] No TAD boundary for " + chr_key + "!")
                n_ad_interaction_ends += 1

        line = fp.readline()
fp.close()
print("\t[INFO] ... done.")

print("\t[INFO] Number of directed interactions ending in directed digests: " + str(n_dd_interaction_ends))
print("\t[INFO] Number of undirected interactions ending in undirected digests: " + str(n_ud_interaction_ends))
print("\t[INFO] Number of directed interactions ending in ambiguous digests: " + str(n_ad_interaction_ends))

print("\t[INFO] Median distance (redundant) between digests involved in directed interactions only and next TAD boundary: " + str(np.median(dd_dist_array_redundant)))
print("\t[INFO] Median distance (redundant) between digests involved in undirected interactions only and next TAD boundary: " + str(np.median(ud_dist_array_redundant)))
print("\t[INFO] Median distance (redundant) between digests involved in directed and undirected interactions and next TAD boundary: " + str(np.median(ad_dist_array_redundant)))

print("\t[INFO] Writing distances to files ...")
tab_stream_dist_to_tad_dd_output = open(out_prefix + "_dist_to_tad_dd_redundant.tsv", 'wt')
for dist in dd_dist_array_redundant:
    tab_stream_dist_to_tad_dd_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_dd_output.close()

tab_stream_dist_to_tad_ud_output = open(out_prefix + "_dist_to_tad_ud_redundant.tsv", 'wt')
for dist in ud_dist_array_redundant:
    tab_stream_dist_to_tad_ud_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_ud_output.close()

tab_stream_dist_to_tad_ad_output = open(out_prefix + "_dist_to_tad_ad_redundant.tsv", 'wt')
for dist in ad_dist_array_redundant:
    tab_stream_dist_to_tad_ad_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_ad_output.close()

print("[INFO] ... done.")


print("[INFO] Writing BED files containing one digest region for each interaction end ...")

bed_stream_dist_to_tad_dd_output = open(out_prefix + "_dist_to_tad_dd_redundant.bed", 'wt')
for digest_region in dd_interaction_end_array_redundant:
    bed_stream_dist_to_tad_dd_output.write(digest_region + '\n')
bed_stream_dist_to_tad_dd_output.close()

bed_stream_dist_to_tad_ud_output = open(out_prefix + "_dist_to_tad_ud_redundant.bed", 'wt')
for digest_region in ud_interaction_end_array_redundant:
    bed_stream_dist_to_tad_ud_output.write(digest_region + '\n')
bed_stream_dist_to_tad_ud_output.close()

bed_stream_dist_to_tad_ad_output = open(out_prefix + "_dist_to_tad_ad_redundant.bed", 'wt')
for digest_region in ad_interaction_end_array_redundant:
    bed_stream_dist_to_tad_ad_output.write(digest_region + '\n')
bed_stream_dist_to_tad_ad_output.close()

print("[INFO] ... done.")

print()
print("Table row:")
print("out_prefix" + '\t' +

      "n_digests_from_dir_inter" + '\t' +       # Number of digests involved in directed interactions
      "n_digests_from_undir_inter" + '\t' +     # Number of digests involved in undirected interactions

      "n_dd_digest_set" + '\t' +                # Number of digests involved in directed interactions only
      "n_ud_digest_set" + '\t' +                # Number of digests involved in undirected interactions only
      "n_ad_digest_set" + '\t' +                # Number of digests involved in directed and undirected interactions

      "median_dist_dd_unique" + '\t' +          # Median distance between directed digests and next TAD
      "median_dist_ud_unique" + '\t' +          # Median distance between undirected digests and next TAD
      "median_dist_ad_unique" + '\t' +          # Median distance between ambiguous digests and next TAD

      "n_dd_interaction_ends" + '\t' +          # Number of directed interactions ending in directed digests
      "n_ud_interaction_ends" + '\t' +          # Number of undirected interactions ending in undirected digests
      "n_ad_interaction_ends" + '\t' +          # Number of undirected interactions ending in ambiguous digests

      "median_dist_dd_redundant" + '\t' +       # Median distance between directed interaction ends and next TAD
      "median_dist_ud_redundant"  + '\t' +      # Median distance between undirected interaction ends and next TAD
      "median_dist_ad_redundant"                # Median distance between ambiguous interaction ends and next TAD
      )
print(out_prefix + '\t' +

      str(len(digests_from_dir_inter)) + '\t' +
      str(len(digests_from_undir_inter)) + '\t' +

      str(len(dd_digest_set)) + '\t' +
      str(len(ud_digest_set)) + '\t' +
      str(len(ad_digest_set)) + '\t' +

      str(np.median(dd_dist_array)) + '\t' +
      str(np.median(ud_dist_array)) + '\t' +
      str(np.median(ad_dist_array)) + '\t' +

      str(n_dd_interaction_ends) + '\t' +
      str(n_ud_interaction_ends) + '\t' +
      str(n_ad_interaction_ends) + '\t' +

      str(np.median(dd_dist_array_redundant)) + '\t' +
      str(np.median(ud_dist_array_redundant)) + '\t' +
      str(np.median(ad_dist_array_redundant))
      )
