#!/usr/bin/env python
import diachrscripts_toolkit as dclass
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean
import time


"""
fdr.py
This script estimates the signficance of the observation of a certain amount of interactions that
are nominally signicant under a binomial test with the null hypothesis that simple and twisted
read pairs are equally probably.
It inputs a capture Hi-C (CHC) dataset (that has been processed by Diachromatic). For each interaction
(one line in the Diachromatic file), it creates a ChcInteraction object and performs a binomial test
with n=total number of interactions, p=0.5, and k=min(# simple read pairs, # twisted read pairs). The result
of this test is assessed as to whether it is significant at a nominal level of alpha=0.05. We then count
the total number of such significant tests (call this number W). Then, we permute the directionality of the data by choosing the 
direction of each read at random and we recalculate the number of "significant" interactions. We do this
n times (n=10000 by default) and record the total number of times that we see W or more "significant interactions 
(call this number w). Our estimated p-value is then w/W.
"""


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine significance of numbers of x:y interactions with x+y=k.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--iter-num', help='Number of iterations.', default=1000)
parser.add_argument('--nominal-alpha', help='Nominal alpha. P-value threshold used for permutation analysis.', default=0.05)
parser.add_argument('--min-digest-dist', help='All interactions with smaller distances will be discarded.', default=20000)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
iter_num = int(args.iter_num)
nominal_alpha = float(args.nominal_alpha)
min_digest_dist = int(args.min_digest_dist)

if diachromatic_interaction_file == None:
    print("--interaction-file option required")
    exit(1)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --iter-num: " + str(iter_num))
print("\t[INFO] --nominal-alpha: " + str(nominal_alpha))
print("\t[INFO] --min-digest-dist: " + str(min_digest_dist))


### Define auxiliary functions
##############################

# Use a dictionary to keep track of P-values.
# key - a string like 2:7
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


def count_significant_pvals_in_permutation(chc_interactions, n_dict, n_alpha=nominal_alpha):
    """
    This function performs the permutation and determines the number of significant permuted interactions for
    one iteration.

    :param chc_interactions: List of original interactions.
    :param n_dict: Dictionary with read pair numbers as keys and interaction numbers as values.
    :param n_alpha: Nominal alpha.
    :return: Number of significant permuted interactions for one iteration.
    """
    n_perm_significant = 0
    random_dict = random_numbers_dict(n_dict)

    for chci in chc_interactions:
        # get permuted counts
        n = chci.n_simple + chci.n_twisted
        n_simple = random_dict[n].pop()
        n_twisted = n - n_simple
        # calculate P-value
        key = "{}-{}".format(n_simple, n_twisted)
        if key in pval_memo:
            pv = pval_memo[key]
        else:
            pv = binomial_p_value(n_simple, n_twisted)
            pval_memo[key] = pv
        # count significant interactions
        if pv <= n_alpha:
            n_perm_significant += 1

    return n_perm_significant


### Start execution
###################

#  Input diachromatic file
n_dict = {} # dictionary that stores the numbers of interactions with n read pairs; is used as input for 'random_numbers_dict()' to generate random count efficiently
chc_interactions = [] # list of original interactions
n_indefinable_cutoff = dclass.find_indefinable_n(nominal_alpha) # minimum number of read pairs required for significance
print("[INFO] Ingesting diachromatic file at", diachromatic_interaction_file, "...")
nsig_o = 0 # number of observed significant interactions
n_interaction = 0
n_cis_long_range_interaction = 0
n_trans_short_range_interaction = 0
n_indefinable_interaction = 0
n_undirected_interaction = 0


with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
    for line in fp:
        n_interaction +=1
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))

        if fields[0] == fields[4] and int(fields[5])-int(fields[2]) >= min_digest_dist: # cis long range

            n_cis_long_range_interaction += 1

            n_simple, n_twisted = fields[8].split(':')
            n_simple = int(n_simple)
            n_twisted = int(n_twisted)
            n = n_simple + n_twisted
            if n < n_indefinable_cutoff: # skip interactions that cannot be significant due to small number of read pairs
                n_indefinable_interaction +=1
                continue

            chci = dclass.Interaction(line)
            chc_interactions.append(chci)

            key = "{}-{}".format(n_simple, n_twisted)
            if key in pval_memo:
                pv = pval_memo[key]
            else:
                pv = chci.get_binomial_p_value()
                pval_memo[key] = pv

            if pv < nominal_alpha:
                nsig_o += 1
            else:
                n_undirected_interaction += 1

            if n in n_dict:
                n_dict[n] +=1
            else:
                n_dict[n] = 1
        else:
            n_trans_short_range_interaction +=1

print("[INFO] Total number interactions: {}".format(n_interaction))
print("[INFO] Total number of (discarded) trans short range interactions: {}".format(n_trans_short_range_interaction))
print("[INFO] Total number of cis long range interactions: {}".format(n_cis_long_range_interaction))
print("[INFO] Number of indefinable cis long range interactions: {}".format(n_indefinable_interaction))
print("[INFO] Number of undirected cis long range interactions: {}".format(n_undirected_interaction))
print("[INFO] Number of cis long range interactions with nominally significant P-values: {}".format(nsig_o))

# Perform permutation anlysis
print("[INFO] Performing permutation analysis with " + str(iter_num) + " iterations ...")
random_better_than_observed = 0
nsig_p_list = [] # stores numbers of significant interactions for each iteration
t = time.process_time()
for n in range(1, iter_num+1):
    nsig_p = count_significant_pvals_in_permutation(chc_interactions, n_dict, nominal_alpha)
    nsig_p_list.append(nsig_p)
    if nsig_p >= nsig_o:
        random_better_than_observed += 1
    if n % int(iter_num / 10) == 0:
        elapsed_time = time.process_time() - t
        print("\t[INFO] " + str(n) + " permutations for " + str(len(chc_interactions)) +  " interactions performed in " + str(elapsed_time) + " sec.")


### Print results
#################

nsig_p_average = mean(nsig_p_list) # calculate average number of significant permuted interactions
percentage_observed = nsig_o/len(chc_interactions)
percentage_permuted = nsig_p_average/len(chc_interactions)

print("OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_TRANS_SHORT_INTERACTION"
               "\tN_CIS_LONG_INTERACTION"
               "\tN_INDEFINABLE_CIS_LONG_INTERACTION"
               "\tN_UNDIRECTED_CIS_LONG_INTERACTION"
               "\tN_DIRECTED_CIS_LONG_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_CIS_LONG_INTERACTION"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED")

print(out_prefix + "\t"
               + str(iter_num) + "\t"
               + str(nominal_alpha) + "\t"
               + str(n_indefinable_cutoff) + "\t"
               + str(n_interaction) + "\t"
               + str(n_trans_short_range_interaction) + "\t"
               + str(n_cis_long_range_interaction) + "\t"
               + str(n_indefinable_interaction) + "\t"
               + str(n_undirected_interaction) + "\t"
               + str(nsig_o) + "\t"
               + str(nsig_p_average) + "\t"
               + str(random_better_than_observed) + "\n")


print("{} out of {} permutations had more signficant p values than in the observed data".format(random_better_than_observed, str(iter_num)))

file_name = out_prefix + "_permutation_summary.txt"
f_output = open(file_name, 'wt')

f_output.write("OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_TRANS_SHORT_INTERACTION"
               "\tN_CIS_LONG_INTERACTION"
               "\tN_INDEFINABLE_CIS_LONG_INTERACTION"
               "\tN_UNDIRECTED_CIS_LONG_INTERACTION"
               "\tN_DIRECTED_CIS_LONG_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_CIS_LONG_INTERACTION"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED\n")

f_output.write(out_prefix + "\t"
               + str(iter_num) + "\t"
               + str(nominal_alpha) + "\t"
               + str(n_indefinable_cutoff) + "\t"
               + str(n_interaction) + "\t"
               + str(n_trans_short_range_interaction) + "\t"
               + str(n_cis_long_range_interaction) + "\t"
               + str(n_indefinable_interaction) + "\t"
               + str(n_undirected_interaction) + "\t"
               + str(nsig_o) + "\t"
               + str(nsig_p_average) + "\t"
               + str(random_better_than_observed) + "\n")


f_output.close()

file_name = out_prefix + "_n_sig_permuted_interactions.txt"
f_output = open(file_name, 'wt')
for nsig_p in nsig_p_list:
    f_output.write(str(nsig_p) + "\n")
f_output.close()
