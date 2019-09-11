#!/usr/bin/env python
import diachrscripts_toolkit as dclass
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean
import time


"""
simulate_k_interactions.py
This script estimates the signficance of the observation of a certain amount of interactions that
are nominally signicant under a binomial test with the null hypothesis that directed and undirected
interactions are equally probably.
It inputs a capture Hi-C (CHC) dataset (that has been processed by Diachromatic). For each interaction
(one line in the Diachromatic file), it creates a ChcInteraction object and performs a binomial test
with n=total number of interactions, p=0.5, and k=observed number of directed interactions. The result
of this test is assessed as to whether it is significant at a nominal level of alpha=0.05. We then count
the total number of such significant tests (call this number W). Then, we permute the directionality of the data by choosing the 
direction of each read at random and we recalculate the number of "significant" interactions. We do this
n times (n=1000 by default) and record the total number of times that we see W or more "significant interactions 
(call this number w). Our estimated p-value is then w/W.
"""


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine significance of numbers of x:y interactions with x+y=k.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--iter-num', help='Number of iterations.', default=10000)
parser.add_argument('--bonferroni-adjustment', help='Perform Bonferroni adjustment of P-value threshold.', default=False, action='store_true')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
ITER_NUM = int(args.iter_num)
bonferroni_adjustment = args.bonferroni_adjustment
diachromatic_interaction_file = args.interaction_file

if diachromatic_interaction_file == None:
    print("--interaction-file option required")
    exit(1)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --iter-num: " + str(ITER_NUM))
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)


# Use a disctionary to keep track of p-values.
# key- a string like 2:7
# value - our corresponding binomial p-value
# note -- use this as a global variable in this script!
pval_memo = defaultdict(float)


class ChcInteraction:
    """
    This class represents one individual interaction, i.e., one line from the diachromatic file
    Note that for the purposes of determining significance, we are ONLY interested in the number
    of simple and twisted read pairs
    """
    def __init__(self, count_string):
        """
        :param count_string: for instance, 0:2, 4:12 etc. n_simple:n_twisted
        """
        fields = count_string.split(':')
        if len(fields) != 2:
            raise TypeError("Malformed counts string", count_string)
        self.n_simple = int(fields[0])
        self.n_twisted = int(fields[1])
        self.N = self.n_simple + self.n_twisted

    def get_counts(self):
        return self.n_simple, self.n_twisted

    def get_permuted_counts(self):
        """
        Choose a number k from [0,N] by random uniform distribution
        :return: k, N-k
        """
        k = binom.rvs(self.N, 0.5, size=1)
        return k, self.N - k

def binomial_p_value(n_simple, n_twisted):
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
    random_numbers_dict = {}
    for n, num in n_dict.items():
        random_numbers_dict[n] = list(binom.rvs(n, p=0.5, size=num))
    return random_numbers_dict


####  Input diachromatic file
#############################

n_dict = {}

chc_interactions = []
observed_significant_interactions = None
permutation_significant_interactions = []
NOMINAL_ALPHA = 0.05
print("[INFO] NOMINAL_ALPHA: " + str(NOMINAL_ALPHA))

# apply Bonferroni correction of P-value threshold
correct_for_multiple_testing = True
if bonferroni_adjustment:

    print("\t[INFO] Adjusting P-value treshold for multiple testing using Bonferroni")

    # determine number of tests to be performed
    n_cis_long_range_interaction = 0
    n_trans_short_range_interaction = 0
    with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
        for line in fp:
            fields = line.rstrip('\n').split('\t')
            if fields[0] == fields[4] and int(fields[5]) - int(fields[2]) > 10000:
                n_cis_long_range_interaction += 1
            else:
                n_trans_short_range_interaction += 1
    fp.close()

    # adjust P-value threshold
    print("\t[INFO] Total number of tests to be performed: " + str(n_cis_long_range_interaction))
    NOMINAL_ALPHA = NOMINAL_ALPHA / n_cis_long_range_interaction
    print("\t[INFO] Adjusted NOMINAL_ALPHA: " + str(NOMINAL_ALPHA))

# iterate interactions
print("[INFO] Ingesting diachromatic file at ", diachromatic_interaction_file, " ...")
nsig_o = 0 # number of observed significant interactions
n_cis_long_range_interaction = 0
n_trans_short_range_interaction = 0
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
    for line in fp:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))

        if fields[0] == fields[4] and int(fields[5])-int(fields[2]) > 10000: # cis long range

            #counts = fields[8]
            chci = dclass.Interaction(line) #ChcInteraction(counts)
            chc_interactions.append(chci)
            n_simple, n_twisted = chci.get_counts()
            key = "{}-{}".format(n_simple, n_twisted)
            if key in pval_memo:
                pv = pval_memo[key]
            else:
                pv = chci.get_binomial_p_value()
                pval_memo[key] = pv

            if pv < NOMINAL_ALPHA:
                nsig_o += 1
            n = n_simple + n_twisted
            if n in n_dict:
                n_dict[n] +=1
            else:
                n_dict[n] = 1

            n_cis_long_range_interaction += 1
        else:
            n_trans_short_range_interaction +=1

print("[INFO] Total number of cis long range interactions: {}".format(n_cis_long_range_interaction  ))
print("[INFO] Total number of trans short range interactions: {}".format(n_trans_short_range_interaction))
print("[INFO] Number of nominally significant P-values: {}".format(nsig_o))
## Now calculate 10,000 permutations and see if we get more p values by chance
def count_signicant_pvals_in_permutation():
    n_perm_significant = 0
    random_dict = random_numbers_dict(n_dict)
    for chci in chc_interactions:
        n = chci.n_simple + chci.n_twisted
        n_simple = random_dict[n].pop()
        n_twisted = n - n_simple
        key = "{}-{}".format(n_simple, n_twisted)
        if key in pval_memo:
            pv = pval_memo[key]
        else:
            pv = chci.get_binomial_p_value()
            pval_memo[key] = pv

        if pv <= NOMINAL_ALPHA:
            n_perm_significant += 1
    return n_perm_significant

random_better_than_observed = 0
nsig_p_list = [] # stores numbers of significant interactions for each iteration
t = time.process_time()
for n in range(ITER_NUM):
    if n % int(ITER_NUM/10) == 0:
        elapsed_time = time.process_time() - t
        print("\t[INFO] " + str(n) + " permuations for " + str(len(chc_interactions)) +  " interactions performed in " + str(elapsed_time) + " sec.")
    nsig_p = count_signicant_pvals_in_permutation()
    nsig_p_list.append(nsig_p)
    permutation_significant_interactions.append(nsig_p)
    if nsig_p >= nsig_o:
        random_better_than_observed += 1

# print numbers to file
nsig_p_average = mean(nsig_p_list)
percentage_observed = nsig_o/len(chc_interactions)
percentage_permuted = nsig_p_average/len(chc_interactions)
print("OUT_PREFIX\tITER_NUM\tTOTAL_INTERACTION_NUM\tNSIG_OBSERVED\tMEAN_NSIG_PERMUATATED\tPERCENTAGE_NSIG_OBSERVED\tPERCENTAGE_MEAN_NSIG_PERMUATATED")
print(out_prefix + "\t" + str(ITER_NUM) + "\t" + str(len(chc_interactions)) + "\t" + str(nsig_o) + "\t" + str(nsig_p_average) + "\t" + str(percentage_observed) + "\t" + str(percentage_permuted))
print("{} out of {} permutations had more signficant p values than in the observed data".format(random_better_than_observed, len(permutation_significant_interactions)))

file_name = out_prefix + "_permuation_analysis_results.txt"
f_output = open(file_name, 'wt')
f_output.write("OUT_PREFIX\tITER_NUM\tNSIG_OBSERVED\tMEAN_NSIG_PERMUATATED\tPERCENTAGE_NSIG_OBSERVED\tPERCENTAGE_MEAN_NSIG_PERMUATATED\n")
f_output.write(out_prefix + "\t" + str(ITER_NUM) + "\t" + str(nsig_o) + "\t" + str(nsig_p_average) + "\t" + str(percentage_observed) + "\t" + str(percentage_permuted))
f_output.close()

file_name = out_prefix + "_nsig_permutated_results.txt"
f_output = open(file_name, 'wt')
for nsig_p in nsig_p_list:
    f_output.write(str(nsig_p) + "\n")
f_output.close()

FDR = nsig_p_average/nsig_o
print(FDR)
