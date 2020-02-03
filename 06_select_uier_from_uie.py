#!/usr/bin/env python

"""
This script takes an enhanced interaction file (created with the script '05_define_di_and_uie.py') with interaction
categories DI, for directed interactions, and UIE, for exclusive undirected interactions (UIE), specified in the third
column and selects reference interactions UIER from UIE that are comparable to DI with respect to enrichment status of
digests, read pair numbers and, optionally, interaction distances.

The script implements a two pass approach:

   1. Iterate interactions and collect information about DI with respect to enrichment status of digests, read pair
   numbers and interaction distances.
   2. Iterate interactions a second time and use the collected information to select reference interactions that are
   comparable to DI.

DI and UIER interactions are written to a file in enhanced interaction format with the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_uier.tsv.gz'

Column 3 of the created file contains either DI or UIER.

Finally, summary statistics ...

"""


import argparse
import gzip
import diachrscripts_toolkit
import numpy


### Parse command line
######################

parser = argparse.ArgumentParser(description='Select reference interactions from exclusive undirected interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file created with \'05_define_di_and_uie.py\'. Column 3 contains either \'DI\' or\'UIE\'.', required=True)
parser.add_argument('--adjust-to-interaction-distance', help='Adjust selected reference interactions to interaction distances of directed interactions.', action='store_true', default=False)
args = parser.parse_args()

out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
adjust_to_interaction_distance = float(args.adjust_to_interaction_distance)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --adjust-to-interaction-distance: " + str(adjust_to_interaction_distance))


### Define auxiliary functions
##############################


### Prepare variables, data structures and streams for output files
###################################################################

# Number of exclusive undirected interactions
dir_inter_num = 0

dir_inter_aa_num = 0
dir_inter_ai_num = 0
dir_inter_ii_num = 0

dir_inter_aa_percentage = None
dir_inter_ai_percentage = None
dir_inter_ii_percentage = None

dir_inter_aa_rp_array = []
dir_inter_ai_rp_array = []
dir_inter_ii_rp_array = []

# Number of exclusive undirected interactions
undir_exc_inter_num = 0

undir_exc_inter_aa_num = 0
undir_exc_inter_ai_num = 0
undir_exc_inter_ii_num = 0

undir_exc_inter_aa_percentage = None
undir_exc_inter_ai_percentage = None
undir_exc_inter_ii_percentage = None

undir_exc_inter_aa_rp_array = []
undir_exc_inter_ai_rp_array = []
undir_exc_inter_ii_rp_array = []

# Prepare stream for output of filtered interactions annotated with respect to exclusive undirected interactions
enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_with_di_and_uier.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')


### 1st pass: Collect information about DI and UIE
##################################################

print("[INFO] 1st pass: Collect information about DI and UIE ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 100000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        if interaction_category == 'DI':
            dir_inter_num += 1
            if enrichment_pair_tag == 'AA':
                dir_inter_aa_num +=1
                dir_inter_aa_rp_array.append(rp_total)
            elif enrichment_pair_tag == 'AI':
                dir_inter_ai_num += 1
                dir_inter_ai_rp_array.append(rp_total)
            elif enrichment_pair_tag == 'II':
                dir_inter_ii_num += 1
                dir_inter_ii_rp_array.append(rp_total)
        elif interaction_category == 'UIE' or interaction_category == 'UII':
            undir_exc_inter_num += 1
            if enrichment_pair_tag == 'AA':
                undir_exc_inter_aa_num +=1
                undir_exc_inter_aa_rp_array.append(rp_total)
            elif enrichment_pair_tag == 'AI':
                undir_exc_inter_ai_num += 1
                undir_exc_inter_ai_rp_array.append(rp_total)
            elif enrichment_pair_tag == 'II':
                undir_exc_inter_ii_num += 1
                undir_exc_inter_ii_rp_array.append(rp_total)
        else:
            print("[Error] Interaction category must be either \'DI\' or \'UIE\'!")
            exit(1)

        line = fp.readline()

dir_inter_aa_percentage = "{0:.2f}".format(100 * dir_inter_aa_num / dir_inter_num)
dir_inter_ai_percentage = "{0:.2f}".format(100 * dir_inter_ai_num / dir_inter_num)
dir_inter_ii_percentage = "{0:.2f}".format(100 * dir_inter_ii_num / dir_inter_num)

# Calculate quantiles and mean
dir_inter_aa_rp_q1 = int(numpy.quantile(dir_inter_aa_rp_array, 0.25))
dir_inter_aa_rp_q2 = int(numpy.quantile(dir_inter_aa_rp_array, 0.5))
dir_inter_aa_rp_q3 = int(numpy.quantile(dir_inter_aa_rp_array, 0.75))
dir_inter_aa_rp_mean = int(numpy.mean(dir_inter_aa_rp_array))

dir_inter_ai_rp_q1 = int(numpy.quantile(dir_inter_ai_rp_array, 0.25))
dir_inter_ai_rp_q2 = int(numpy.quantile(dir_inter_ai_rp_array, 0.5))
dir_inter_ai_rp_q3 = int(numpy.quantile(dir_inter_ai_rp_array, 0.75))
dir_inter_ai_rp_mean = int(numpy.mean(dir_inter_ai_rp_array))

dir_inter_ii_rp_q1 = int(numpy.quantile(dir_inter_ii_rp_array, 0.25))
dir_inter_ii_rp_q2 = int(numpy.quantile(dir_inter_ii_rp_array, 0.5))
dir_inter_ii_rp_q3 = int(numpy.quantile(dir_inter_ii_rp_array, 0.75))
dir_inter_ii_rp_mean = int(numpy.mean(dir_inter_ii_rp_array))

undir_exc_inter_aa_percentage = "{0:.2f}".format(100 * undir_exc_inter_aa_num / undir_exc_inter_num)
undir_exc_inter_ai_percentage = "{0:.2f}".format(100 * undir_exc_inter_ai_num / undir_exc_inter_num)
undir_exc_inter_ii_percentage = "{0:.2f}".format(100 * undir_exc_inter_ii_num / undir_exc_inter_num)

# Calculate quantiles and mean
undir_exc_inter_aa_rp_q1 = int(numpy.quantile(undir_exc_inter_aa_rp_array, 0.25))
undir_exc_inter_aa_rp_q2 = int(numpy.quantile(undir_exc_inter_aa_rp_array, 0.5))
undir_exc_inter_aa_rp_q3 = int(numpy.quantile(undir_exc_inter_aa_rp_array, 0.75))
undir_exc_inter_aa_rp_mean = int(numpy.mean(undir_exc_inter_aa_rp_array))

undir_exc_inter_ai_rp_q1 = int(numpy.quantile(undir_exc_inter_ai_rp_array, 0.25))
undir_exc_inter_ai_rp_q2 = int(numpy.quantile(undir_exc_inter_ai_rp_array, 0.5))
undir_exc_inter_ai_rp_q3 = int(numpy.quantile(undir_exc_inter_ai_rp_array, 0.75))
undir_exc_inter_ai_rp_mean = int(numpy.mean(undir_exc_inter_ai_rp_array))

undir_exc_inter_ii_rp_q1 = int(numpy.quantile(undir_exc_inter_ii_rp_array, 0.25))
undir_exc_inter_ii_rp_q2 = int(numpy.quantile(undir_exc_inter_ii_rp_array, 0.5))
undir_exc_inter_ii_rp_q3 = int(numpy.quantile(undir_exc_inter_ii_rp_array, 0.75))
undir_exc_inter_ii_rp_mean = int(numpy.mean(undir_exc_inter_ii_rp_array))


print()
print("Total number of directed interactions: " + str(dir_inter_num))
print("\tWithin 'AA': " + str(dir_inter_aa_num) + " (" + str(dir_inter_aa_percentage) + "%)")
print("\tWithin 'AI': " + str(dir_inter_ai_num) + " (" + str(dir_inter_ai_percentage) + "%)")
print("\tWithin 'II': " + str(dir_inter_ii_num) + " (" + str(dir_inter_ii_percentage) + "%)")
print()
print("\tAA: Q1 = " + str(dir_inter_aa_rp_q1) + ", " + "Q2 = " + str(dir_inter_aa_rp_q2) + ", " + "Q3 = " + str(dir_inter_aa_rp_q3) + ", " + "Mean = " + str(dir_inter_aa_rp_mean))
print("\tAI: Q1 = " + str(dir_inter_ai_rp_q1) + ", " + "Q2 = " + str(dir_inter_ai_rp_q2) + ", " + "Q3 = " + str(dir_inter_ai_rp_q3) + ", " + "Mean = " + str(dir_inter_ai_rp_mean))
print("\tII: Q1 = " + str(dir_inter_ii_rp_q1) + ", " + "Q2 = " + str(dir_inter_ii_rp_q2) + ", " + "Q3 = " + str(dir_inter_ii_rp_q3) + ", " + "Mean = " + str(dir_inter_ii_rp_mean))



print()
print("Number of exclusive undirected interactions: " + str(undir_exc_inter_num))
print("\tWithin 'AA': " + str(undir_exc_inter_aa_num) + " (" + str(undir_exc_inter_aa_percentage) + "%)")
print("\tWithin 'AI': " + str(undir_exc_inter_ai_num) + " (" + str(undir_exc_inter_ai_percentage) + "%)")
print("\tWithin 'II': " + str(undir_exc_inter_ii_num) + " (" + str(undir_exc_inter_ii_percentage) + "%)")
print()
print("\tAA: Q1 = " + str(undir_exc_inter_aa_rp_q1) + ", " + "Q2 = " + str(undir_exc_inter_aa_rp_q2) + ", " + "Q3 = " + str(undir_exc_inter_aa_rp_q3) + ", " + "Mean = " + str(undir_exc_inter_aa_rp_mean))
print("\tAI: Q1 = " + str(undir_exc_inter_ai_rp_q1) + ", " + "Q2 = " + str(undir_exc_inter_ai_rp_q2) + ", " + "Q3 = " + str(undir_exc_inter_ai_rp_q3) + ", " + "Mean = " + str(undir_exc_inter_ai_rp_mean))
print("\tII: Q1 = " + str(undir_exc_inter_ii_rp_q1) + ", " + "Q2 = " + str(undir_exc_inter_ii_rp_q2) + ", " + "Q3 = " + str(undir_exc_inter_ii_rp_q3) + ", " + "Mean = " + str(undir_exc_inter_ii_rp_mean))


