#!/usr/bin/env python

"""
This script takes a path to a Diachromatic interaction file and determines the overall significance of directed
interactions  by randomization of simple and twisted read pairs.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean, stdev
from numpy import exp, log
import multiprocessing as mp
import itertools
import numpy as np

from diachr.random_permutation import RandomPermutation
from diachr.binomial_model import BinomialModel


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine overall significance of directed interactions by randomization of simple and twisted read pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i','--interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('-n','--iter-num', help='Number of iterations.', default=1000)
parser.add_argument('-a','--nominal-alpha', help='Nominal alpha. P-value threshold used to define significant interactions.', default=0.05)
parser.add_argument('-s','--random-seed', help='Seed for randomization of simple and twisted read pairs.', default=42)
parser.add_argument('-t','--thread-num', help='Number of threads.', default=1)

parser.add_argument('-m','--usemod', help="Use new module", dest='usemod', action='store_true')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
ITER_NUM = int(args.iter_num)
NOMINAL_ALPHA = float(args.nominal_alpha)
NLN_NOMINAL_ALPHA = -log(NOMINAL_ALPHA)
random_seed = int(args.random_seed)
thread_num = int(args.thread_num)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --iter-num: " + str(ITER_NUM))
print("\t[INFO] --nominal-alpha: " + str(NOMINAL_ALPHA))
print("\t\t[INFO] Negative of natural logarithm: " + str(NLN_NOMINAL_ALPHA))
print("\t[INFO] --random-seed: " + str(random_seed))
print("\t[INFO] --thread-num: " + str(thread_num))

if args.usemod:
    print("[INFO] Using module")
    randomP = RandomPermutation(combined_interaction_file=diachromatic_interaction_file, alpha=NOMINAL_ALPHA, threads=thread_num,
                    prefix=out_prefix)
    randomP.permute(ITER_NUM=ITER_NUM)
    exit(1)


### Define auxiliary functions
##############################

# Dictionary that keeps track of already calculated P-values
#    key - a string like 2-7
#    value - our corresponding binomial p-value
#    note -- use this as a global variable in this script!
pval_memo = defaultdict(float)

p_values = BinomialModel()
# Unchanged
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0

# Using P-value class with exp(-nln_pval)/2
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	316119	40203	19646.36	127.98	160.62	0

# Using P-value class with exp(-nln_pval)
# results/tmp_results/02/02	100	0.05	5	1000000	>643678<	332283	24039	8053.22	91.56	174.60	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	332283	24039	8053.22	91.56	174.60	0
# results/tmp_results/02/02	100	0.05	5	1000000	643678	332283	24039	8053.22	91.56	174.60	0

# Using P-value class with exp(-nln_pval) -> adjusted indef-function to two-sided test
# results/tmp_results/02/02	100	0.05	6	1000000	706490	269471	24039	8034.44	75.56	211.80	0

# Using nln(p_val_thresh) instead of p_val_thresh -> Should not change results
# results/tmp_results/02/02	100	0.05	6	1000000	706490	269471	24039	8034.44	75.56	211.80	0
# results/tmp_results/02/02	100	0.05	6	1000000	706490	269471	24039	8034.44	75.56	211.80	0
# results/tmp_results/02/02	100	0.05	6	1000000	706490	269471	24039	8034.44	75.56	211.80	0

# Set random seed
np.random.seed(random_seed)

def perform_one_iteration():
    """
    This function performs one iteration of the randomization analysis. It randomizes the simple and twisted read
    pair counts of all interactions and then determines and returns the number of significant interactions.

    :return: Number of significant interactions after randomization
    """

    # Init number of significant interactions after randomization
    n_random_significant = 0

    # Iterate dictionary with numbers of interactions foreach read pair number n
    for n, i_num in N_DICT.items():

        # Generate random simple read pair counts for current n
        simple_count_list = list(binom.rvs(n, p = 0.5, size = i_num))

        for simple_count in simple_count_list:

            # Get twisted count
            twisted_count = n - simple_count

            # Get binomial P-value
            nln_pv = p_values.get_binomial_nnl_p_value(simple_count, twisted_count)

            # Count significant interactions
            if NLN_NOMINAL_ALPHA < nln_pv:
                n_random_significant += 1

    return n_random_significant


def perform_m_iterations(m_iter):

    w_sig_p_list = []

    print("\t[INFO] Batch: Performing " + str(m_iter) + " iterations ...")

    for i in range(1, m_iter + 1):
        w_sig_p_list.append(perform_one_iteration())

    return w_sig_p_list


### Prepare variables, data structures and streams for output files
###################################################################

# Dictionary that stores the numbers of interactions with n read pairs
N_DICT = {}

# Minimum number of read pairs required for significance
smallest_sig_n, pv_indefinable_cutoff = p_values.find_smallest_significant_n(NOMINAL_ALPHA)

# Total number of input interactions
n_interaction = 0

# Number of interaction that cannot be significant
n_indefinable_interaction = 0

# Number of significant interactions observed for original data (W)
n_directed_interaction = 0

# Number of not significant interactions observed for original data
n_undirected_interaction = 0

# Prepare stream for output of summary statistics
tab_file_stats_output = out_prefix + "_permutation_stats.txt"
tab_stream_stats_output = open(tab_file_stats_output, 'wt')

# Prepare stream for output of significant interaction numbers for all iterations (w_i)
txt_file_w_sig_p_output = out_prefix + "_permutation_w_sig_p.txt"
txt_stream_w_sig_p_output = open(txt_file_w_sig_p_output, 'wt')


### Start execution
###################

print("[INFO] Iterating Diachromatic interaction file ...")
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    for line in fp:

        # Count total number of interactions
        n_interaction += 1

        # Report progress
        if n_interaction % 1000000 == 0:
            print("\t\t[INFO]", n_interaction, "interactions processed ...")

        # Split Diachromatic interaction line
        fields = line.rstrip('\n').split('\t')

        # Check format of Diachromatic interaction line
        if len(fields) < 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))

        # Extract simple and twisted read pair counts from Diachromatic interaction line
        n_simple, n_twisted = fields[8].split(':')
        n_simple = int(n_simple)
        n_twisted = int(n_twisted)
        n = n_simple + n_twisted

        # Skip interactions that cannot be significant
        if n < smallest_sig_n:
            n_indefinable_interaction += 1
            continue

        # Get P-value of interaction
        nln_pv = p_values.get_binomial_nnl_p_value(n_simple, n_twisted)


        # Count interaction as directed or undirected
        if NLN_NOMINAL_ALPHA < nln_pv:
            n_directed_interaction += 1
        else:
            n_undirected_interaction += 1

        # Add the sum of simple and twisted read pair counts to dictionary that will be used for randomization
        if n in N_DICT:
            N_DICT[n] +=1
        else:
            N_DICT[n] = 1


# Output some statistics about original interactions
print("[INFO] Total number interactions: {}".format(n_interaction))
print("\t[INFO] Number of indefinable interactions (discarded): {}".format(n_indefinable_interaction))
print("\t[INFO] Number of undirected interactions: {}".format(n_undirected_interaction))
print("\t[INFO] Number of interactions with nominally significant P-values (<{}): {}".format(NOMINAL_ALPHA, n_directed_interaction))

# Perform permutation analysis
print("[INFO] Performing permutation analysis with " + str(ITER_NUM) + " iterations ...")

# Init pool with 'thread_num' processes
pool = mp.Pool(processes = thread_num)

# Perform permutation for 'thread_num' batches with balanced numbers of iterations
batch_iter_nums = list(itertools.repeat(int(ITER_NUM / thread_num), thread_num))
if 0 < ITER_NUM - thread_num * int(ITER_NUM / thread_num):
    # If there is a remainder, add it to the first element of list
    batch_iter_nums[0] = batch_iter_nums[0] + ITER_NUM - thread_num * int(ITER_NUM / thread_num)
results = [pool.apply_async(perform_m_iterations, args=(batch_iter_num,)) for batch_iter_num in batch_iter_nums]

# Combine results from different batches
batch_results = [p.get() for p in results]
print("[INFO] Combining results from different batches ...")

# List with numbers of significant interactions after permutation for each iteration
w_sig_p_list = list(itertools.chain.from_iterable(batch_results))

# Get number of iterations with more significant interactions as compared to the original data
random_better_than_observed = len([i for i in w_sig_p_list if i > n_directed_interaction])

# Calculate average number of significant permuted interactions
w_sig_p_average = mean(w_sig_p_list)

# Calculate standard deviation of significant permuted interactions
w_sig_p_sd = stdev(w_sig_p_list)

# Compute Z-score
z_score = "{0:.2f}".format((n_directed_interaction - w_sig_p_average) / w_sig_p_sd)

# Calculate percentage of significant interactions for original data
percentage_observed = "{0:.2f}".format(n_directed_interaction / (n_directed_interaction + n_undirected_interaction))

# Calculate percentage of significant interactions for the permuted data
percentage_permuted = "{0:.2f}".format(w_sig_p_average / (n_directed_interaction + n_undirected_interaction))


### Print summary statistics to screen
######################################

print("\t[INFO] OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_INDEFINABLE_INTERACTION"
               "\tN_UNDIRECTED_INTERACTION"
               "\tN_DIRECTED_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_INTERACTION"
               "\tSD_PERMUTED_DIRECTED_INTERACTION"
               "\tZ_SCORE"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED")

print("\t[INFO] " + out_prefix + "\t"
      + str(ITER_NUM) + "\t"
      + str(NOMINAL_ALPHA) + "\t"
      + str(smallest_sig_n) + "\t"
      + str(n_interaction) + "\t"
      + str(n_indefinable_interaction) + "\t"
      + str(n_undirected_interaction) + "\t"
      + str(n_directed_interaction) + "\t"
      + "{0:.2f}".format(w_sig_p_average) + "\t"
      + "{0:.2f}".format(w_sig_p_sd) + "\t"
      + str(z_score) + "\t"
      + str(random_better_than_observed) + "\n")

print("{} out of {} permutations had more signficant p values than in the observed data.".format(random_better_than_observed, str(ITER_NUM)))


### Write results to the file
#############################

# Write table with parameters and summary statistics to file
tab_stream_stats_output.write("OUT_PREFIX"
               "\tITER_NUM"
               "\tNOMINAL_ALPHA"
               "\tINDEF_RP_CUTOFF"               
               "\tN_INTERACTION"
               "\tN_INDEFINABLE_INTERACTION"
               "\tN_UNDIRECTED_INTERACTION"
               "\tN_DIRECTED_INTERACTION"
               "\tMEAN_PERMUTED_DIRECTED_INTERACTION"
               "\tSD_PERMUTED_DIRECTED_INTERACTION"
               "\tZ_SCORE"
               "\tN_PERMUTED_BETTER_THAN_OBSERVED\n")

tab_stream_stats_output.write(out_prefix + "\t"
                              + str(ITER_NUM) + "\t"
                              + str(NOMINAL_ALPHA) + "\t"
                              + str(smallest_sig_n) + "\t"
                              + str(n_interaction) + "\t"
                              + str(n_indefinable_interaction) + "\t"
                              + str(n_undirected_interaction) + "\t"
                              + str(n_directed_interaction) + "\t"
                              + "{0:.2f}".format(w_sig_p_average) + "\t"
                              + "{0:.2f}".format(w_sig_p_sd) + "\t"
                              + str(z_score) + "\t"
                              + str(random_better_than_observed) + "\n")

tab_stream_stats_output.close()

# Write numbers of significant interactions for all iterations to file
for w_sig_p in w_sig_p_list:
    txt_stream_w_sig_p_output.write(str(w_sig_p) + "\n")
txt_stream_w_sig_p_output.close()
