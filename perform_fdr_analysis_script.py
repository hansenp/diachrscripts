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
perform_fdr_analysis_script.py
This script uses a simple procedure to estimate the FDR of directed interactions.

"""

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine a P-value threshold that corresponds to a given FDR threshold.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--fdr-threshold', help='Use this switch to estimate a P-value cutoff that corresponds to a given FDR threshold.', default=0.25)
parser.add_argument('--p-val-c-min', help='Smallest P-value cutoff.', default=0.00025)
parser.add_argument('--p-val-c-max', help='Largest P-value cutoff.', default=0.05)
parser.add_argument('--p-val-step-size', help='P-value step size.', default=0.00025)
parser.add_argument('--min-digest-dist', help='All interactions with smaller distances will be discarded.', default=20000)


args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
fdr_threshold = float(args.fdr_threshold)
p_val_c_min = float(args.p_val_c_min)
p_val_c_max = float(args.p_val_c_max)
p_val_step_size = float(args.p_val_step_size)
min_digest_dist = int(args.min_digest_dist)

if diachromatic_interaction_file == None:
    print("--interaction-file option required")
    exit(1)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --fdr-threshold: " + str(fdr_threshold))
print("\t[INFO] --p-val-c-min: " + str(p_val_c_min))
print("\t[INFO] --p-val-c-max: " + str(p_val_c_max))
print("\t[INFO] --p-val-step-size: " + str(p_val_step_size))
print("\t[INFO] --min-digest-dist: " + str(min_digest_dist))


### Define auxiliary functions
##############################

# Use a disctionary to keep track of p-values.
# key- a string like 2:7
# value - our corresponding binomial p-value
# note -- use this as a global variable in this script!
pval_memo = defaultdict(float)
def binomial_p_value(n_simple, n_twisted):
    """
    Locally defined method for the calculation of the binomial P-value that uses a dictionary that keeps track of
    P-values that already have been calculated.

    :param n_simple: number of simple read pairs
    :param n_twisted: Number of twisted read pairs
    :return: Binomial P-value
    """
    key = "{}-{}".format(n_simple, n_twisted)
    if key in pval_memo:
        return pval_memo[key]
    else:
        if n_simple < n_twisted:
            p_value = 1 - binom.cdf(n_twisted - 1, n_simple + n_twisted, 0.5)
            pval_memo[key] = p_value
            return p_value
        else:
            p_value = 1 - binom.cdf(n_simple - 1, n_simple + n_twisted, 0.5)
            pval_memo[key] = p_value
            return p_value


def random_numbers_dict(n_dict):
    """
    This function can be used to generate all random numbers required for one permutation at once,
    which is far more efficient that generating random numbers one by one.

    :param n_dict: Dictionary with read pair numbers as keys and interaction numbers as values.
    :return: Dictionary with read pair numbers as keys and lists of random numbers of simple read pairs as values.
    """
    random_numbers_dict = {}
    for n, num in n_dict.items():
        random_numbers_dict[n] = list(binom.rvs(n, p=0.5, size=num))
    return random_numbers_dict


def get_pvals_permuted_counts(chc_interactions, n_dict):
    """
    This function is used for the FDR analysis only.

    :param n_dict: Dictionary with read pair numbers as keys and interaction numbers as values.
    :return: List of P-values for permuted interactions.
    """
    pvals_permuted_counts = []
    random_dict = random_numbers_dict(n_dict)
    for chci in chc_interactions:
        n = chci.n_simple + chci.n_twisted
        n_simple = random_dict[n].pop()
        n_twisted = n - n_simple
        pvals_permuted_counts.append(binomial_p_value(n_simple, n_twisted))
    return pvals_permuted_counts


### Start execution
###################

#  Input diachromatic file
print("[INFO] Ingesting diachromatic file at ", diachromatic_interaction_file, " ...")
n_dict = {} # dictionary that stores the numbers of interactions with n read pairs; is used as input for 'random_numbers_dict()' to generate random count efficiently
chc_interactions = [] # list of original interactions
n_cis_long_range_interaction = 0
n_trans_short_range_interaction = 0
p_val_o_list = [] # list of P-values for observed interactions
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
    for line in fp:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))

        if fields[0] == fields[4] and int(fields[5])-int(fields[2]) >= min_digest_dist: # cis long range

            chci = dclass.Interaction(line)
            chc_interactions.append(chci)
            n_simple, n_twisted = chci.get_counts()
            key = "{}-{}".format(n_simple, n_twisted)
            if key in pval_memo:
                pv = pval_memo[key]
            else:
                pv = chci.get_binomial_p_value()
                pval_memo[key] = pv

            p_val_o_list.append(pv)
            n = n_simple + n_twisted

            if n in n_dict:
                n_dict[n] +=1
            else:
                n_dict[n] = 1

            n_cis_long_range_interaction += 1
        else:
            print(line)
            n_trans_short_range_interaction +=1

print("[INFO] Total number of cis long range interactions: {}".format(n_cis_long_range_interaction))
print("[INFO] Total number of trans short range interactions: {}".format(n_trans_short_range_interaction))


# Perform FDR anlysis
print("[INFO] Determining P-value threshold ...")
file_name = out_prefix + "_fdr_analysis_results.txt"
f_output = open(file_name, 'wt')
f_output.write("OUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O" + "\n")
p_val_p_list = get_pvals_permuted_counts(chc_interactions, n_dict)
for pc in np.arange(p_val_c_min, p_val_c_max, p_val_step_size):
    nsig_o = (p_val_o_list <= pc).sum()
    nsig_p = (p_val_p_list <= pc).sum()
    fdr = nsig_p / nsig_o
    f_output.write(out_prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(nsig_p) + "\t" + str(nsig_o) + "\n")
    print("\t" + out_prefix + "\t" + str(fdr) + "\t" + str(pc) + "\t" + str(nsig_p) + "\t" + str(nsig_o))
    if fdr < fdr_threshold:
        fdr_last = fdr
        nsig_o_last = nsig_o
        nsig_p_last = nsig_p
        pc_last = pc
f_output.close()

print("\tOUT_PREFIX\tFDR\tPC\tNSIG_P\tNSIG_O")
print(
    "\t" + out_prefix + "\t" + str(fdr_last) + "\t" + str(pc_last) + "\t" + str(nsig_p_last) + "\t" + str(nsig_o_last))
