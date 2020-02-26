#!/usr/bin/env python
import diachrscripts_toolkit as dclass
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean
import time
import numpy as np


"""
Add P-values and category (S, T, U, UR, NA) to interactions.
"""

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine a P-value threshold that corresponds to a given FDR threshold.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--p-value-cutoff', help='Use this P-value for categorization.', default=0.05)
parser.add_argument('--min-digest-dist', help='All interactions with smaller distances will be discarded.', default=20000)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
p_value_cutoff = float(args.p_value_cutoff)
min_digest_dist = int(args.min_digest_dist)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + diachromatic_interaction_file)
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))
print("\t[INFO] --min-digest-dist: " + str(min_digest_dist))


### Start execution
###################

file_name_original = out_prefix + "_annotated_interactions.tsv"
f_output_original = open(file_name_original, 'wt')

n_interaction_total = 0
n_trans_short_range_interaction = 0
n_non_promoter_promoter_interaction = 0
n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_undirected_reference_interaction = 0
n_indefinable_interaction = 0

# Given the P-value cutoff, find the smallest n that yields a significant P-value
n_indefinable_cutoff = dclass.find_indefinable_n(p_value_cutoff)

# Determine distribution of n for directed interactions
n_dict = dclass.get_n_dict(diachromatic_interaction_file, "ALL", min_digest_dist, p_value_cutoff)

pval_dict = {}

# iterate interactions
print("[INFO] Determining pair category for each interaction in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, 'rt') as fp:

    line = fp.readline().rstrip()

    while line:

        n_interaction_total += 1
        if n_interaction_total%1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        fields = line.split("\t")

        if fields[0] != fields[4]:
            continue # cis

        if int(fields[5]) - int(fields[2]) < 20000:
            continue # short range

        # get simple and twisted read pair counts
        fields2 = fields[8].split(":")
        n_simple = int(fields2[0])
        n_twisted = int(fields2[1])
        n_total = n_simple + n_twisted

        if fields[8] in pval_dict:
            pval = pval_dict[fields[8]]
        else:
            pval = dclass.calculate_binomial_p_value(n_simple, n_twisted)
            pval_dict[fields[8]] = float(pval)

        if n_simple + n_twisted < n_indefinable_cutoff:
            f_output_original.write(line + "\t" + "NA" + "\t" + str(pval) + "\n")
        else:
            if pval < p_value_cutoff:
                if n_twisted < n_simple:
                    f_output_original.write(line + "\t" + "S" + "\t" + str(pval) + "\n")
                else:
                    f_output_original.write(line + "\t" + "T" + "\t" + str(pval) + "\n")
            else:
                if n_total in n_dict and 0 < n_dict[n_total]:
                    f_output_original.write(line + "\t" + "UR" + "\t" + str(pval) + "\n")
                    n_dict[n_total] = n_dict[n_total] - 1
                else:
                    f_output_original.write(line + "\t" + "U" + "\t" + str(pval) + "\n")

        line = fp.readline().rstrip()

    # check whether some reference interactions are missing
    for x in n_dict:
        if  0 < n_dict[x]:
            print("For " + str(n_dict[x]) + " directed interactions with " + str(x) + " read pairs no undirected reference interaction could be slected.")

fp.close()
f_output_original.close()

