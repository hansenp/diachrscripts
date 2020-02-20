#!/usr/bin/env python

"""
This script estimates the significance of the observation of a certain amount of interactions that are nominally
significant under a binomial test with the null hypothesis that simple and twisted read pairs are equally probably.
It inputs a capture Hi-C (CHC) dataset (that has been processed by Diachromatic).

For each interaction (one line in the Diachromatic file), it performs a binomial test with n = total number of
interactions, p = 0.5 and k = min(# simple read pairs (s), # twisted read pairs (t)). The result of this test is
assessed as to whether it is significant at a nominal level of alpha = 0.05. We then count the total number of such
significant tests (call this number W).

Then, we randomize the input by permuting the counts of simple and twisted read pairs of each interaction using
Binom(n = s + t, p = 0.5) and recalculate the number of 'significant' interactions. We do this m times
(m = 1000 by default) and record the number of significant interactions for each iteration (w_i). If we do not observe
an iteration with W < w_i, we estimate the empirical P-value with p < 1/m.

Finally, in order to compare the overall significance of for different datasets, we compute the Z-score from the mean
and standard deviation of w_i.
"""


import diachrscripts_toolkit as dclass
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean, stdev
import multiprocessing as mp
import itertools


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine overall significance of directed interactions using random'
                                             'permutation of simple and twisted read pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--iter-num', help='Number of iterations.', default=1000)
parser.add_argument('--nominal-alpha', help='Nominal alpha. P-value threshold used for permutation analysis.', default=0.05)
parser.add_argument('--thread-num', help='Number of threads.', default=1)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
ITER_NUM = int(args.iter_num)
NOMINAL_ALPHA = float(args.nominal_alpha)
thread_num = int(args.thread_num)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --iter-num: " + str(ITER_NUM))
print("\t[INFO] --nominal-alpha: " + str(NOMINAL_ALPHA))
print("\t[INFO] --thread-num: " + str(thread_num))


### Define auxiliary functions
##############################

# Dictionary that keeps track of already calculated P-values
#    key - a string like 2-7
#    value - our corresponding binomial p-value
#    note -- use this as a global variable in this script!
pval_memo = defaultdict(float)

def binomial_p_value(simple_count, twisted_count):
    """
    Locally defined method for the calculation of the binomial P-value that uses a dictionary that keeps track of
    P-values that already have been calculated.

    :param simple_count: Number of simple read pairs
    :param twisted_count: Number of twisted read pairs
    :return: Binomial P-value
    """

    # Create key from simple and twisted read pair counts
    key = "{}-{}".format(simple_count, twisted_count)

    # Check whether a P-value for this combination of simple and twisted counts has been calculated already
    if key in pval_memo:
        return pval_memo[key]
    else:

        # Calculate P-value and add to dictionary
        if simple_count < twisted_count:
            p_value = 1 - binom.cdf(twisted_count - 1, simple_count + twisted_count, 0.5)
            pval_memo[key] = p_value
            return p_value
        else:
            p_value = 1 - binom.cdf(simple_count - 1, simple_count + twisted_count, 0.5)
            pval_memo[key] = p_value
            return p_value


def perform_one_iteration():
    """
    This function performs one iteration of the random permutation analysis. It permutes the simple and twisted read
    pair counts of all interactions and then determines and returns the number of significant interactions.

    :return: Number of significant interactions after permutation
    """

    # Init number of significant interactions after permutation
    w_perm_significant = 0

    # Iterate dictionary with numbers of interactions foreach read pair number n
    for n, i_num in N_DICT.items():

        # Generate random simple read pair counts for current n
        simple_count_list = list(binom.rvs(n, p = 0.5, size = i_num))

        for simple_count in simple_count_list:

            # Get twisted count
            twisted_count = n - simple_count

            # Get binomial P-value
            key = "{}-{}".format(simple_count, twisted_count)
            if key in pval_memo:
                pv = pval_memo[key]
            else:
                pv = binomial_p_value(simple_count, twisted_count)
                pval_memo[key] = pv

            # Count significant interactions
            if pv < NOMINAL_ALPHA:
                w_perm_significant += 1

    return w_perm_significant


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
n_indefinable_cutoff = dclass.find_indefinable_n(NOMINAL_ALPHA)

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
        if n < n_indefinable_cutoff:
            n_indefinable_interaction += 1
            continue

        # Get P-value of interaction
        key = "{}-{}".format(n_simple, n_twisted)
        if key in pval_memo:
            pv = pval_memo[key]
        else:
            pv = binomial_p_value(n_simple, n_twisted)
            pval_memo[key] = pv

        # Count interaction as directed or undirected
        if pv < NOMINAL_ALPHA:
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

print("OUT_PREFIX"
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

print(out_prefix + "\t"
      + str(ITER_NUM) + "\t"
      + str(NOMINAL_ALPHA) + "\t"
      + str(n_indefinable_cutoff) + "\t"
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
                              + str(n_indefinable_cutoff) + "\t"
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
