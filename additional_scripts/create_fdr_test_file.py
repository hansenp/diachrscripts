#!/usr/bin/env python

"""
With this script the test file for the FDR procedure was generated.
"""

import argparse
from numpy import log, arange
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

parser = argparse.ArgumentParser(description='Evaluate and categorize interactions and select unndirected reference interactions.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('--p-value-max', help='Upper limit for P-value ranges.', default=0.05)
parser.add_argument('--p-value-step', help='Size of P-value ranges.', default=0.00025)
parser.add_argument('--i-count-per-range', help='Size of P-value ranges.', default=10)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.diachromatic_interaction_file
p_value_max = float(args.p_value_max)
p_value_step = float(args.p_value_step)
i_count_per_range = int(args.i_count_per_range)

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file + '\n'
parameter_info += "\t[INFO] --p-value-max: " + str(p_value_max) + '\n'
parameter_info += "\t[INFO] --p-value-step: " + str(p_value_step) + '\n'
parameter_info += "\t[INFO] --i-count-per-range: " + str(i_count_per_range) + '\n'
print(parameter_info)

# Open stream for output interaction file
out_fname = out_prefix + '_diachromatic_fdr_test_file.tsv'
out_fh = open(out_fname, 'wt')

# Load interactions
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()
print(read_file_info_report)

# Evaluate interactions
interaction_set.evaluate_and_categorize_interactions((p_value_max + 2*p_value_step), verbose=True)
eval_cat_info_report = interaction_set.get_eval_cat_info_report()
print(eval_cat_info_report)

# Select interactions in P-value ranges
i_count = 0
p_threshs = arange(p_value_step, p_value_max + p_value_step, p_value_step)
print(p_threshs)
for p_thresh in p_threshs:
    i_count_range = 0
    for d_inter in interaction_set.interaction_list:
        if (p_thresh - p_value_step < d_inter.get_pval()) and d_inter.get_pval() <= p_thresh:
            i_count_range  += 1
            out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
            i_count += 1
        if i_count_range == i_count_per_range:
            break
    if i_count_range < i_count_per_range:
        print("[WARNING] Could not select the required number (only "
              + str(i_count_range) + " of " + str(i_count_per_range) +
              ") of interactions for the P-value range ]"
              + str(p_thresh - p_value_step) + ';' + str(str(p_thresh)) + ']')
    print(str(p_thresh - p_value_step) + '\t' + str(p_thresh) + '\t' + str(i_count_range) + '\t' + str(i_count))

# Open stream for output interaction file
out_fh.close()
print("[INFO] Wrote " + str(i_count) + " interactions to: " + out_fname)
