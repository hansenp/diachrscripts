#!/usr/bin/env python
import argparse
import gzip
from scipy.stats import binom
from collections import defaultdict
import random


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
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file

if  diachromatic_interaction_file == None:
    print("--interaction-file option required")
    exit(1)

print("[INFO] --out_prefix: " + out_prefix)
print("[INFO] --interaction-file: " + diachromatic_interaction_file)


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
        #k = random.randint(0, self.N)
        k = binom.rvs(self.N, 0.5, size=1)
        #return self.n_simple, self.n_twisted
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
# iterate interactions
print("[INFO] Ingesting diachromatic file at ", diachromatic_interaction_file, " ...")
n_significant = 0
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
    for line in fp:
        #print(line)
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            raise TypeError("Malformed diachromatic input line {} (number of fields {})".format(line, len(fields)))
        counts = fields[8]
        chci = ChcInteraction(counts)
        chc_interactions.append(chci)
        n_simple, n_twisted = chci.get_counts()
        pv = binomial_p_value(n_simple, n_twisted)
        if pv < NOMINAL_ALPHA:
            n_significant += 1
        n = n_simple + n_twisted
        if n in n_dict:
            n_dict[n] +=1
        else:
            n_dict[n] = 1
print("[INFO]] Number of nominally significant pvalues: {}".format(n_significant))


## Now calculate 10,000 permutations and see if we get more p values by chance
def count_signicant_pvals_in_permutation():
    n_perm_significant = 0
    random_dict = random_numbers_dict(n_dict)
    for chci in chc_interactions:
        n = chci.n_simple + chci.n_twisted
        n_simple = random_dict[n].pop()
        n_twisted = n - n_simple
        #n_simple, n_twisted = chci.get_permuted_counts()
        pv = binomial_p_value(n_simple, n_twisted)
        if pv < NOMINAL_ALPHA:
            n_perm_significant += 1
    return n_perm_significant


random_better_than_observed = 0
for n in range(100):
    nsig = count_signicant_pvals_in_permutation()
    permutation_significant_interactions.append(nsig)
    print("nsig (permuted): {} nsig (observed: {})".format(nsig, n_significant))
    if nsig >= n_significant:
        random_better_than_observed += 1

print("{} out of {} permutations had more signficant p values than in the observed data".format(random_better_than_observed, len(permutation_significant_interactions)))
