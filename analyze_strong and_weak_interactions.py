#!/usr/bin/env python
from collections import defaultdict

import argparse
from scipy.stats import uniform, binom
import numpy as np
import matplotlib.pyplot as plt

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine significance of numbers of x:y interactions with x+y=k.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--i-num', help='Number of simulated interactions')
parser.add_argument('--n-max', help='Simulate interactions with 1 to n read pairs.')
parser.add_argument('--p-value-cutoff', help='Nominal alpha. P-value threshold used for the classification of directed and undirected interactions.', default=0.05)

args = parser.parse_args()
out_prefix = args.out_prefix
n_max = int(args.n_max)
i_num = int(args.i_num)
p_value_cutoff = float(args.p_value_cutoff)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --i-num: " + str(i_num))
print("\t[INFO] --n-max: " + str(n_max))
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))


### Define auxiliary functions
##############################

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


### Start execution
###################

print("[INFO] " + "Generating numbers of simple and twisted read pairs ...")
n_dict = {} # dictionary that stores the numbers of interactions with n read pairs; is used as input for 'random_numbers_dict()' to generate random counts efficiently
random_n_vec = np.random.randint(low = 1, high = n_max+1, size=i_num) # n has a uniform distribution
for n in random_n_vec:
    if n in n_dict:
        n_dict[n] += 1
    else:
        n_dict[n] = 1
random_numbers_dict = random_numbers_dict(n_dict)
# 'random_numbers_dict' now contains arrays of random simple read numbers for each n

for N in random_numbers_dict:
    #print("N: " + str(len(random_numbers_dict[N]))+ "\t" + str(N) + "\t" + str(random_numbers_dict[N]))
    arr = [0] * (N+1)
    for n_simple in random_numbers_dict[N]:
        arr[n_simple] += 1

    print(str(N) + ": " + str(arr))

print("[INFO] " + "Counting significant interactions for each n ...")
signum_list = [0] * (n_max + 1)
for N in range(1, n_max + 1):
    print(N)
    if N in random_numbers_dict:
        for n_simple in random_numbers_dict[N]:
            n_twisted = N - n_simple
            p_val = binomial_p_value(n_simple, n_twisted)
            if p_val <= p_value_cutoff:
                signum_list[N] += 1

plt.plot(signum_list)
plt.grid(True)
plt.xlabel("N read pairs")
plt.ylabel("#Significant interactions")
sub_title = "#Interactions: " + str(i_num) + " | Max n: " + str(n_max) + " | P-value cutoff: " + str(p_value_cutoff)
plt.suptitle(sub_title)
figure_name = out_prefix + "_uniform_n_sig_interaction_dist.pdf"
plt.savefig(figure_name, format = "pdf")
