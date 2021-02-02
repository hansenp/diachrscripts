#!/usr/bin/env python

"""
This script takes a path to a Diachromatic interaction file and determines the overall significance of directed
interactions  by randomization of simple and twisted read pairs.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
from scipy.stats import binom
from collections import defaultdict
from statistics import mean, stdev
from numpy import log
import multiprocessing as mp
import itertools
import numpy as np
import sys

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet


# Parse command line
####################

parser = argparse.ArgumentParser(description='Determine overall significance of directed interactions by randomization of simple and twisted read pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i','--interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('-n','--iter-num', help='Number of iterations.', default=1000)
parser.add_argument('-a','--nominal-alpha', help='Nominal alpha. P-value threshold used to define significant interactions.', default=0.05)
parser.add_argument('-s','--random-seed', help='Seed for randomization of simple and twisted read pairs.', default=None)
parser.add_argument('-t','--thread-num', help='Number of threads.', default=1)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
iter_num = int(args.iter_num)
nominal_alpha = float(args.nominal_alpha)
nnl_nominal_alpha = -log(nominal_alpha)
if args.random_seed is None:
    random_seed = args.random_seed
else:
    random_seed = int(args.random_seed)
thread_num = int(args.thread_num)

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --interaction-file: " + diachromatic_interaction_file + '\n'
parameter_info += "\t[INFO] --iter-num: " + str(iter_num) + '\n'
parameter_info += "\t[INFO] --nominal-alpha: " + str(nominal_alpha) + '\n'
parameter_info += "\t[INFO] --random-seed: " + str(random_seed) + '\n'
parameter_info += "\t[INFO] --thread-num: " + str(thread_num) + '\n'

print(parameter_info)


# Perform analysis
##################

# Load interactions
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()

# Create RandomizeInteractionSet object and perform analysis
randomize_interactions = RandomizeInteractionSet(interaction_set=interaction_set, random_seed=random_seed)
randomize_interactions_info_dict = randomize_interactions.perform_randomization_analysis(
    nominal_alpha = nominal_alpha,
    iter_num = iter_num,
    thread_num = thread_num,
    verbose = True)


# Create plot and file with summary statistics
##############################################

# Create plot with a histogram for the randomized significant interaction numbers from all iterations
randomize_interactions.get_randomization_info_plot(
    pdf_file_name = out_prefix + '_randomization_analysis.pdf',
    analysis_name = out_prefix)

f_name_summary = out_prefix + "_randomization_analysis_summary.txt"
out_fh = open(f_name_summary, 'wt')

# Chosen parameters
out_fh.write(parameter_info + '\n')

# Report on reading files
out_fh.write(randomize_interactions.get_randomization_info_report() + '\n')

out_fh.close()

sys.exit(0)

### Define auxiliary functions
##############################

# Dictionary that keeps track of already calculated P-values
#    key - a string like 2-7
#    value - our corresponding binomial p-value
#    note -- use this as a global variable in this script!
pval_memo = defaultdict(float)

p_values = BinomialModel()

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
    for n, i_num in RP_I_DICT.items():

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

# Dictionary that stores for each total read pair number (n) the number of interactions
RP_I_DICT = {}

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

# Get list of Diachromatic interaction objects
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file)
d_inter_list = interaction_set.interaction_list

# Iterate Dichromatic interaction objects
print("[INFO] Iterating list Diachromatic interaction objects ...")
for d_inter in d_inter_list:

    # Count total number of interactions
    n_interaction += 1

    # Report progress
    if n_interaction % 1000000 == 0:
        print("\t\t[INFO]", n_interaction, "interactions processed ...")

    # Skip interactions that cannot be significant
    if d_inter.rp_total < smallest_sig_n:
        n_indefinable_interaction += 1
        continue

    # Get P-value of interaction
    nln_pv = p_values.get_binomial_nnl_p_value(d_inter.n_simple, d_inter.n_twisted)

    # Count directed and undirected interactions
    if NLN_NOMINAL_ALPHA < nln_pv:
        n_directed_interaction += 1
    else:
        n_undirected_interaction += 1

    # Add total number of read pair counts to dictionary that will be used for randomization
    if d_inter.rp_total in RP_I_DICT:
        RP_I_DICT[d_inter.rp_total] += 1
    else:
        RP_I_DICT[d_inter.rp_total] = 1

print("[INFO] ... done.")


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
