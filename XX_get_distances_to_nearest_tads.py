#!/usr/bin/env python

"""
This script takes as input an enhanced interaction file and a BED file with TAD regions. Two digest sets are derived
from the interaction file:

   1. Digests that are involved in at least one directed interaction.
   2. Digest that are involved in undirected reference interactions and in no directed interaction.

For each digest, the distance to the nearest TAD boundary is determined and added to alist that corresponds to one of
the two digest sets.

The TAD regions are stored in an object of the class 'TADs' that is essentially a hash map. The names of the
chromosomes serve as keys and for each key there is a sorted list of boundary coordinates. The class has a function that
returns the coordinate of the next TAD boundary for a given coordinate. This is accomplished through binary search.
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

args = parser.parse_args()

out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
tad_region_bed_file = args.tad_region_bed_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --tad-region-bed-file: " + tad_region_bed_file)


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
            #print(chr_key + "\t" + str(self.chr_dict[chr_key]))


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

print("[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        x = line.split('\t')
        n_simple = x[4].split(':')[0]
        n_twisted = x[4].split(':')[1]

        # Add digest to set for directed interactions
        if interaction_category == "DIII" or interaction_category == "DIII" or interaction_category == "DIII":
            digests_from_dir_inter.add(chr_a + '\t' + str(sta_a) + '\t' + str(end_a))

        # Add digest to set for undirected interactions
        if interaction_category == "UIRII" or interaction_category == "UIRII" or interaction_category == "UIRII":
            digests_from_undir_inter.add(chr_b + '\t' + str(sta_b) + '\t' + str(end_b))

        line = fp.readline()

fp.close()
print("[INFO] ... done.")

digests_from_undir_inter_wo_digests_from_dir_inter = digests_from_undir_inter - digests_from_dir_inter
print("Number of digests from directed interactions: " + str(len(digests_from_dir_inter)))
print("Number of digests from undirected interactions: " + str(len(digests_from_undir_inter_wo_digests_from_dir_inter)))


total_dist_dir = 0
dir_dist_array = []
for digest in digests_from_dir_inter:
    chr_key, sta, end = digest.split('\t')
    digest_center_pos = int(sta) + int((int(end) - int(sta))/2)
    dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_key,digest_center_pos)
    if dist == -1:
        print("Warning: No TAD boundary for " + chr_key + "!")
    else:
        total_dist_dir = total_dist_dir + dist
        dir_dist_array.append(dist)

print("Average distance of digests from directed interactions: " + str(total_dist_dir/len(digests_from_dir_inter)))


total_dist_undir = 0
undir_dist_array = []
for digest in digests_from_undir_inter_wo_digests_from_dir_inter:
    chr_key, sta, end = digest.split('\t')
    digest_center_pos = int(sta) + int((int(end) - int(sta))/2)
    dist = TAD_boundaries.get_distance_to_nearest_tad_boundary(chr_key,digest_center_pos)
    if dist == -1:
        print("Warning: No TAD boundary for " + chr_key + "!")
    else:
        total_dist_undir = total_dist_undir + dist
        undir_dist_array.append(dist)

print("Average distance of digests from undirected interactions: " + str(total_dist_undir/len(digests_from_undir_inter_wo_digests_from_dir_inter)))

print(np.median(dir_dist_array))
print(np.median(undir_dist_array))


# Write distances to files

tab_stream_dist_to_tad_dir_output = open(out_prefix + "_dist_to_tad_dir.tsv", 'wt')
for dist in dir_dist_array:
    tab_stream_dist_to_tad_dir_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_dir_output.close()

tab_stream_dist_to_tad_undir_output = open(out_prefix + "_dist_to_tad_undir.tsv", 'wt')
for dist in undir_dist_array:
    tab_stream_dist_to_tad_undir_output.write(str(dist) + '\n')
tab_stream_dist_to_tad_undir_output.close()


