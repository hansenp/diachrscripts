#!/usr/bin/env python
import argparse
import gzip
import diachrscripts_toolkit as dclass
from decimal import *
from math import isinf

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine numbers of x:y interactions with x+y=k.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file

print("[INFO] --out_prefix: " + out_prefix)
print("[INFO] --interaction-file: " + diachromatic_interaction_file)


### Define auxiliary classes and functions
##########################################

class kInteractionCounter:

    def __init__(self, max_k=1000):
        self.k_interaction_dict = {}
        for i in range(2, max_k + 1):
            self.max_k = max_k
            self.k_interaction_dict[i] = [0, 0, 0] # [NK, NK_ZERO_SIMPLE, NK_ZERO_TWISTED]
            # NK - all interactions with k read pairs
            # NK_ZERO_SIMPLE - interactions with k read pairs and zero simple read pairs
            # NK_ZERO_TWISTED - interactions with k read pairs and zero twisted read pairs

    def increment_k_interaction(self, n_simple_rp, n_twisted_rp):
        if (n_simple_rp + n_twisted_rp) < self.max_k:
            self.k_interaction_dict[n_simple_rp + n_twisted_rp][0] += 1
            if n_simple_rp == 0:
                self.k_interaction_dict[n_simple_rp + n_twisted_rp][1] += 1
            if n_twisted_rp == 0:
                self.k_interaction_dict[n_simple_rp + n_twisted_rp][2] += 1
        else:
            self.k_interaction_dict[self.max_k][0] += 1
            if n_simple_rp == 0:
                self.k_interaction_dict[self.max_k][1] += 1
            if n_twisted_rp == 0:
                self.k_interaction_dict[self.max_k][2] += 1

    def get_fraction_zero_interactions(self, k):
        n_zero_interactions = self.k_interaction_dict[k][1] + self.k_interaction_dict[k][2]
        if self.k_interaction_dict[k][0] == 0:
            return str(0.0)
        else:
            return str(Decimal(float(1.0*n_zero_interactions/self.k_interaction_dict[k][0])).quantize(Decimal(10) ** -2))

    def get_expected_number_zero_interactions(self, k):
        p = 2 * (0.5 ** k)
        return p * (self.k_interaction_dict[k][0])

    def print_k_interaction_counts_to_file(self, file_path_name):

        # open file handle
        f = open(file_path_name, "w+")

        # print header to file
        f.write("K\tN_INTER\tN_ZERO_SIMP\tN_ZERO_TWIST\tN_ZERO_TOT\tF_ZERO_TOT\tN_ZERO_TOT_EXP\tF_EXP_AMONG_OBS\tLOG10_P_VAL\n")
        for i in range(2, self.max_k + 1):

            # extract and format output
            k = str(i)
            n_interaction = str(self.k_interaction_dict[i][0])
            n_zero_simple = str(self.k_interaction_dict[i][1])
            n_zero_twisted = str(self.k_interaction_dict[i][2])
            n_zero_total = str(self.k_interaction_dict[i][1] + self.k_interaction_dict[i][2])
            f_zero_total = dclass.get_string_formatted_fraction(self.k_interaction_dict[i][1] + self.k_interaction_dict[i][2], self.k_interaction_dict[i][0])
            n_zero_total_expected = str(int(round(self.get_expected_number_zero_interactions(i))))
            frac_expected_among_observed = dclass.get_string_formatted_fraction(self.get_expected_number_zero_interactions(i), self.k_interaction_dict[i][1] + self.k_interaction_dict[i][2])
            p_value, error_code = dclass.get_x_inter_binomial_log10_p_value(self.k_interaction_dict[i][1] + self.k_interaction_dict[i][2], self.k_interaction_dict[i][0], i)

            if error_code == 1:
                p_value == "NA"
            elif error_code == 2:
                p_value = "<" + str(round(p_value,2))
            else:
                p_value = str(round(p_value,2))

            # print to output to file
            f.write(k + "\t" + n_interaction\
                    + "\t" + n_zero_simple\
                    + "\t" + n_zero_twisted\
                    + "\t" + n_zero_total\
                    + "\t" + f_zero_total\
                    + "\t" + n_zero_total_expected \
                    + "\t" + frac_expected_among_observed \
                    + "\t" + p_value \
                    + "\n")

        # close file handle
        f.close()


### Start execution
###################

k_interaction_counter = kInteractionCounter(100)

n_interaction_total = 0
n_trans_short_range_interaction = 0

# iterate interactions
print("[INFO] Determining number of k-interactions with zero simple or twisted read pairs in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    line = fp.readline()

    while line:

        if n_interaction_total%1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(10000)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        # increment interaction counts
        k_interaction_counter.increment_k_interaction(interaction.n_simple, interaction.n_twisted)

        line = fp.readline()

    print("... done.")

fp.close()

print("[INFO] Total number of interactions: " + str(n_interaction_total))
print("[INFO] Number of trans and short range interactions: " + str(n_trans_short_range_interaction) + " (discarded)")

file_path_name = out_prefix + "_k_inter.tab"
k_interaction_counter.print_k_interaction_counts_to_file(file_path_name = file_path_name)
print("[INFO] Output written to: " + file_path_name)
print("[INFO] Done.")
