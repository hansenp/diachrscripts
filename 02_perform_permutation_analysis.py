#!/usr/bin/env python

"""
This script estimates the signficance of the observation of a certain amount of interactions that are nominally
significant under a binomial test with the null hypothesis that simple and twisted read pairs are equally probably.
It inputs a capture Hi-C (CHC) dataset (that has been processed by Diachromatic).

For each interaction (one line in the Diachromatic file), it creates a ChcInteraction object and performs a binomial
test with n = total number of interactions, p = 0.5 and k = min(# simple read pairs (s), # twisted read pairs (t)).
The result of this test is assessed as to whether it is significant at a nominal level of alpha = 0.05. We then count
the total number of such significant tests (call this number W).

Then, we randomize the input by permuting the counts of simple and twisted read pairs of each interaction using
Binom(n = s + t, p = 0.5) and recalculate the number of 'significant' interactions.

We do this m times (m = 1000 by default) and record the total number of times that we see W or more significant
interactions (call this number w). Our estimated p-value is then w/W.
"""


import diachrscripts_toolkit as dclass
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
from statistics import mean
import time
import multiprocessing as mp


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

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --iter-num: " + str(iter_num))
print("\t[INFO] --nominal-alpha: " + str(nominal_alpha))
print("\t[INFO] --min-digest-dist: " + str(min_digest_dist))


### Define auxiliary functions
##############################

# Dictionary that keeps track of already calculated P-values
#    key - a string like 2-7
#    value - our corresponding binomial p-value
#    note -- use this as a global variable in this script!
pval_memo = defaultdict(float)

def binomial_p_value(n_simple, n_twisted):
    """
    Locally defined method for the calculation of the binomial P-value that uses a dictionary that keeps track of
    P-values that already have been calculated.

    :param n_simple: number of simple read pairs
    :param n_twisted: Number of twisted read pairs
    :return: Binomial P-value
    """

    # Create key from simple and twisted read pair counts
    key = "{}-{}".format(n_simple, n_twisted)

    # Check whether a P-value for this combination of simple and twisted counts has been calculated already
    if key in pval_memo:
        return pval_memo[key]
    else:

        # Calculate P-value and add to dictionary
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

    :param n_dict: Dictionary with read pair numbers as keys and interaction numbers as values
    :return: Dictionary with read pair numbers as keys and lists of random numbers of simple read pairs as values
    """

    # Init dictionary that will be returned
    random_numbers_dict = {}

    # Iterate input dictionary with interaction numbers for each n
    for n, num in n_dict.items():
        # Generate a list of random simple read pair counts for all interactions with n read pairs
        random_numbers_dict[n] = list(binom.rvs(n, p=0.5, size=num))

    # Return dictionary containing simulated numbers for all interactions
    return random_numbers_dict


def count_significant_pvals_in_permutation(chc_interactions, n_dict, n_alpha=nominal_alpha):
    """
    This function performs one iteration of the random permutation analysis. It permutes the simple and twisted read
    pair counts of all interactions and then determines the number of significant interactions.

    :param chc_interactions: List of original interactions
    :param n_dict: Dictionary with read pair numbers as keys and interaction numbers as values
    :param n_alpha: Nominal alpha
    :return: Number of significant interactions after permutation
    """

    # Init number of significant interactions after permutation
    n_perm_significant = 0

    # Generate random numbers for all interactions at once
    random_dict = random_numbers_dict(n_dict)

    # Iterate interaction objects
    for chci in chc_interactions:

        # Get permuted counts of simple and twisted read pairs
        n = chci.n_simple + chci.n_twisted
        n_simple = random_dict[n].pop()
        n_twisted = n - n_simple

        # Get binomial P-value
        key = "{}-{}".format(n_simple, n_twisted)
        if key in pval_memo:
            pv = pval_memo[key]
        else:
            pv = binomial_p_value(n_simple, n_twisted)
            pval_memo[key] = pv

        # Count significant interactions
        if pv <= n_alpha:
            n_perm_significant += 1

    return n_perm_significant


def perform_n_iterations_of_permuatations(n_iter, chc_interactions, n_dict, n_alpha):

    n_sig_p_list = []

    for n in range(1, n_iter + 1):
        n_sig_p_list.append(count_significant_pvals_in_permutation(chc_interactions, n_dict, n_alpha))

    return n_sig_p_list


### Prepare variables, data structures and streams for output files
###################################################################

# Dictionary that stores the numbers of interactions with n read pairs
n_dict = {}

# List of objects for original interactions
chc_interactions = []

# Minimum number of read pairs required for significance
n_indefinable_cutoff = dclass.find_indefinable_n(nominal_alpha)

# Total number of input interactions
n_interaction = 0

# Total number of cis long range input interactions
n_cis_long_range_interaction = 0

# Total number of trans short range input interactions
n_trans_short_range_interaction = 0

# Number of interaction that cannot be significant
n_indefinable_interaction = 0

# Number of significant interactions observed for original data
n_directed_interaction = 0

# Number of not significant interactions observed for original data
n_undirected_interaction = 0


### Start execution
###################

print("[INFO] Iterating Diachromatic interaction file ...")
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    for line in fp:

        # Count total number of interactions
        n_interaction +=1

        if n_interaction % 100000 == 0:
            print("\t\t[INFO]", n_interaction, "interactions processed ...")

        # Split Diachromatic interaction line
        fields = line.rstrip('\n').split('\t')

        # Check format of Diachromatic interaction line
        if len(fields) < 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))

        # Filter for cis long range interactions
        if fields[0] == fields[4] and int(fields[5])-int(fields[2]) >= min_digest_dist: # cis long range

            # Count cis long range interactions
            n_cis_long_range_interaction += 1

            # Extract simple and twisted read pair counts from Diachromatic interaction line
            n_simple, n_twisted = fields[8].split(':')
            n_simple = int(n_simple)
            n_twisted = int(n_twisted)
            n = n_simple + n_twisted

            # Skip interactions that cannot be significant
            if n < n_indefinable_cutoff:
                n_indefinable_interaction +=1
                continue

            # Create interaction object and append list of interactions
            chci = dclass.Interaction(line)
            chc_interactions.append(chci)

            # Get P-value of interaction
            key = "{}-{}".format(n_simple, n_twisted)
            if key in pval_memo:
                pv = pval_memo[key]
            else:
                pv = chci.get_binomial_p_value()
                pval_memo[key] = pv

            # Count interaction as directed or undirected
            if pv < nominal_alpha:
                n_directed_interaction += 1
            else:
                n_undirected_interaction += 1

            # Add the sum of simple and twisted read pair counts to dictionary required for efficient randomization
            if n in n_dict:
                n_dict[n] +=1
            else:
                n_dict[n] = 1

        else:
            # Count trans short range interactions
            n_trans_short_range_interaction +=1


# Output some statistics about original interactions
print("[INFO] Total number interactions: {}".format(n_interaction))
print("[INFO] Total number of (discarded) trans short range interactions: {}".format(n_trans_short_range_interaction))
print("[INFO] Total number of cis long range interactions: {}".format(n_cis_long_range_interaction))
print("\t[INFO] Number of indefinable interactions: {}".format(n_indefinable_interaction))
print("\t[INFO] Number of undirected interactions: {}".format(n_undirected_interaction))
print("\t[INFO] Number of interactions with nominally significant P-values: {}".format(n_directed_interaction))

# Perform permutation analysis
print("[INFO] Performing permutation analysis with " + str(iter_num) + " iterations ...")

# Number of iterations with more significant interactions as compared to the original data
random_better_than_observed = 0

# List with numbers of significant interactions after permutation for each iteration
nsig_p_list = []

# Record start time
t = time.process_time()

# Init pool with 2 processes
pool = mp.Pool(processes=4)

n_threads = 8

results = [pool.apply(perform_n_iterations_of_permuatations, args=(10, chc_interactions, n_dict, nominal_alpha)) for x in range(1,10 + 1)]
print(len(results))
print(len(results[0]))
elapsed_time = time.process_time() - t
print(elapsed_time)

# Perform random permutations
for n in range(1, iter_num + 1):

    # Permute counts of all interactions and determine the number of significant interactions
    nsig_p = count_significant_pvals_in_permutation(chc_interactions, n_dict, nominal_alpha)

    # Append number of significant interactions for this iteration to list
    nsig_p_list.append(nsig_p)

    # Count the number of iterations with more significant interactions as compared to the original data
    if nsig_p >= n_directed_interaction:
        random_better_than_observed += 1

    # Report progress
    if n % int(iter_num / 10) == 0:
        elapsed_time = time.process_time() - t
        print("\t[INFO] " + str(n) + " permutations for " + str(len(chc_interactions)) +  " interactions performed in " + str(elapsed_time) + " sec.")

# Calculate average number of significant permuted interactions
nsig_p_average = mean(nsig_p_list)

# Calculate percentage of significant interactions for original data
percentage_observed = n_directed_interaction / len(chc_interactions)

# Calculate percentage of significant interactions for the permuted data
percentage_permuted = nsig_p_average/len(chc_interactions)


### Print results to screen
###########################

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
      + str(n_directed_interaction) + "\t"
      + str(nsig_p_average) + "\t"
      + str(random_better_than_observed) + "\n")


print("{} out of {} permutations had more signficant p values than in the observed data".format(random_better_than_observed, str(iter_num)))


### Write results to the file
#############################

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
               + str(n_directed_interaction) + "\t"
               + str(nsig_p_average) + "\t"
               + str(random_better_than_observed) + "\n")


f_output.close()

file_name = out_prefix + "_n_sig_permuted_interactions.txt"
f_output = open(file_name, 'wt')
for nsig_p in nsig_p_list:
    f_output.write(str(nsig_p) + "\n")
f_output.close()
