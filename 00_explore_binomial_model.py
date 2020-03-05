#!/usr/bin/env python

"""
This script is to investigate the relationship between signal strength, measured as the total number of read pairs,
and the binomial model for directionality of interactions.

For this purpose, a specified number of interactions is simulated, whereby the the total numbers of read pairs are drawn
from a uniform distribution. For individual interactions with n read pairs, the number of simple read pairs is drawn
from a binomial distribution with p = 0.5 (null model) and the number of twisted read pairs is set to n minus the number
of simple read pairs. Subsequently, the simulated interactions are evaluated for statistical significance given the null
model and a specified significance threshold. Finally, a plot for the number of significant simulated interactions for
each n is generated and, in addition, the corresponding numbers are written to a tab separated text file.

If a Diachromatic interaction files is specified, the numbers of significant empirical interactions for each n will also
be determined and written to file.

The two tab separated output text files can be used to reproduce Figure S1 in R.
"""

from collections import defaultdict
import argparse
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Explore binomial model using simuated and real interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--i-num', help='Number of simulated interactions')
parser.add_argument('--n-max', help='Simulate interactions with 1 to n read pairs.')
parser.add_argument('--p-value-cutoff', help='P-value threshold used for the classification of directed and undirected interactions.', default=0.05)
parser.add_argument('--diachromatic-interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
n_max = int(args.n_max)
i_num = int(args.i_num)
p_value_cutoff = float(args.p_value_cutoff)
diachromatic_interaction_file = args.diachromatic_interaction_file


print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --i-num: " + str(i_num))
print("\t[INFO] --n-max: " + str(n_max))
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))
if diachromatic_interaction_file != None:
    print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)


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
        p_value = dclass.calculate_binomial_p_value(simple_count, twisted_count)
        pval_memo[key] = p_value
        return p_value


### Prepare variables, data structures and streams for output files
###################################################################

# Dictionary that stores the numbers of significant interactions for each n
N_DICT_SIM = {}
N_DICT_EMP = {}

# List containing counts of significant simple and twisted interactions for each n
signum_list = [0] * (n_max + 1)

# Name for PDF file with significant simulated interactions versus n plot
pdf_name = out_prefix + "_sig_interactions_vs_uniform_n.pdf"

# Name for text file with significant empirical interactions for each n
sim_tab_file_name = out_prefix + "_sig_interactions_vs_uniform_n.tab"
sim_tab_stream_name = open(sim_tab_file_name, 'wt')

if diachromatic_interaction_file != None:

    # Name for text file with significant empirical interactions for each n
    emp_tab_file_name = out_prefix + "_sig_interactions_vs_empirical_n.tab"
    emp_tab_stream_name = open(emp_tab_file_name, 'wt')


### Start execution
###################

print("[INFO] " + "Generating random numbers of simple and twisted read pairs ...")

# Random vector for n with uniform distribution
random_n_vec = np.random.randint(low = 1, high = n_max  + 1, size = i_num)
for n in random_n_vec:
    if n in N_DICT_SIM:
        N_DICT_SIM[n] += 1
    else:
        N_DICT_SIM[n] = 1


print("[INFO] " + "Counting significant interactions for each n ...")

# Iterate dictionary with numbers of interactions for each read pair number n
for n, i in N_DICT_SIM.items():

    # Generate random simple read pair counts for current n
    simple_count_list = list(binom.rvs(n, p = 0.5, size = i))

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

        # Count significant interactions for current n
        if pv <= p_value_cutoff:
            signum_list[n] += 1


print("[INFO] " + "Generating plot: Significant simulated interactions versus n ...")
plt.plot(signum_list)
plt.grid(True)
plt.xlabel("n")
plt.ylabel("# Significant interactions")
sub_title = "# Interactions: " + str(i_num) + " | Max n: " + str(n_max) + " | P-value cutoff: " + str(p_value_cutoff)
plt.suptitle(sub_title)
plt.savefig(pdf_name, format = "pdf")


print("[INFO] " + "Writing numbers of significant simulated interactions for each n to text file ...")
for n in range(0, n_max + 1):
    try:
        sim_tab_stream_name.write(str(n) + "\t" + str(signum_list[n]) + "\t" + str(N_DICT_SIM[n]) + "\n")
    except KeyError:
        sim_tab_stream_name.write(str(n) + "\t" + str(signum_list[n]) + "\t" + str(0) + "\n")
sim_tab_stream_name.close()


if diachromatic_interaction_file != None:

    # Count significant interactions in empirical data for each n
    N_DICT_EMP = dclass.get_n_dict(diachromatic_interaction_file, 'ALL', 20000, 0.0039)

    print("[INFO] " + "Writing numbers of significant empirical interactions for each n to text file ...")
    for n in range(0, 2000):
        try:
            emp_tab_stream_name.write(str(n) + "\t" + str(N_DICT_EMP[n]) + "\n")
        except KeyError:
            emp_tab_stream_name.write(str(n) + "\t" + str(0) + "\n")

    emp_tab_stream_name.close()
