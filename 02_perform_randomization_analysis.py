#!/usr/bin/env python

"""
This script takes a path to a Diachromatic interaction file and determines the overall significance of directed
interactions by randomization of simple and twisted read pairs.

The individual steps that are carried out in this script are demonstrated in the following Jupyter notebook:

       diachrscripts/jupyter_notebooks/perform_randomization_analysis.ipynb
"""

import argparse
import sys

from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet


# Parse command line
####################

parser = argparse.ArgumentParser(description='Determine overall significance of directed interactions by '
                                             'randomization of simple and twisted read pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUT_PREFIX')
parser.add_argument('-i','--interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('-n','--iter-num', help='Number of iterations.', default=1000)
parser.add_argument('-a','--nominal-alpha', help='Nominal alpha. P-value threshold used to define significant '
                                                 'interactions.', default=0.05)
parser.add_argument('-x','--nominal-alpha-max', help='Maximal nominal alpha used to define significant interactions.',
                    default=0.05)
parser.add_argument('-y','--nominal-alpha-step', help='Increment by which nominal alphas are increased.',
                    default=0.0025)
parser.add_argument('-s','--random-seed', help='Seed for randomization of simple and twisted read pairs.', default=None)
parser.add_argument('-t','--thread-num', help='Number of threads.', default=0)
parser.add_argument('--analysis-name', help='Name for the analysis shown in the plot.', default='ANALYSIS_NAME')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
iter_num = int(args.iter_num)
nominal_alpha = float(args.nominal_alpha)
nominal_alpha_max = float(args.nominal_alpha_max)
nominal_alpha_step = float(args.nominal_alpha_step)
if args.random_seed is None:
    random_seed = args.random_seed
else:
    random_seed = int(args.random_seed)
thread_num = int(args.thread_num)
analysis_name = args.analysis_name

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --interaction-file: " + diachromatic_interaction_file + '\n'
parameter_info += "\t[INFO] --iter-num: " + str(iter_num) + '\n'
parameter_info += "\t[INFO] --nominal-alpha: " + str(nominal_alpha) + '\n'
parameter_info += "\t[INFO] --random-seed: " + str(random_seed) + '\n'
parameter_info += "\t[INFO] --thread-num: " + str(thread_num) + '\n'
parameter_info += "\t[INFO] --analysis-name: " + analysis_name + '\n'

print(parameter_info)


# Perform analysis
##################

# Load interactions
interaction_set = DiachromaticInteractionSet()
min_rp_num, min_rp_num_pval = interaction_set._p_values.find_smallest_significant_n(nominal_alpha)
interaction_set.parse_file(diachromatic_interaction_file, min_rp_num=min_rp_num, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()
print()

# Create RandomizeInteractionSet object and perform analysis
randomize_interactions = RandomizeInteractionSet(random_seed=random_seed)
randomize_interactions_info_dict = randomize_interactions.perform_randomization_analysis(
    interaction_set=interaction_set,
    nominal_alphas = [nominal_alpha],
    iter_num = iter_num,
    thread_num = thread_num,
    verbose = True)
print()

# Create plot and file with summary statistics
##############################################

# Create plot with a histogram for the randomized significant interaction numbers from all iterations
f_name_plot = out_prefix + '_randomization_analysis.pdf'
randomize_interactions.get_randomization_info_plot(
    pdf_file_name = f_name_plot,
    analysis_name = analysis_name)

f_name_summary = out_prefix + "_randomization_analysis_summary.txt"
out_fh = open(f_name_summary, 'wt')

# Chosen parameters
out_fh.write(parameter_info + '\n')

# Report on reading files
out_fh.write(interaction_set.get_read_file_info_report() + '\n')

# Report on randomization analysis
out_fh.write(randomize_interactions.get_randomization_info_report() + '\n')
out_fh.write(randomize_interactions.get_randomization_info_table_row(out_prefix)+ '\n')

# Report on generated files
generated_file_info = "[INFO] Generated files:" + '\n'
generated_file_info += "\t[INFO] " + f_name_plot + '\n'
generated_file_info += "\t[INFO] " + f_name_summary + '\n'
out_fh.write(generated_file_info + '\n')

out_fh.close()

print(generated_file_info)
