#!/usr/bin/env python

"""
This script takes a file in Diachromatic interaction format and an FDR or P-value threshold, calculates the P-values of
interactions, assigns interactions to directed (DI) and undirected interactions (UI), selects undirected reference
interactions (UIR) from UI and creates a new file with two additional columns for P-value and interaction category.

You can find documentation on this script in the relevant section in the RTD of this repository.

The individual steps that are carried out in this script are demonstrated in the following Jupyter notebook:

       diachrscripts/jupyter_notebooks/evaluate_and_categorize_interactions.ipynb
"""

import argparse
from numpy import arange, log
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet
from diachr.randomize_interaction_set import RandomizeInteractionSet

### Parse command line
######################

parser = argparse.ArgumentParser(description='Evaluate and categorize interactions and select undirected reference '
                                             'interactions.')
parser.add_argument('-o', '--out-prefix', help='Prefix for output.', default='OUT_PREFIX')
parser.add_argument('-d', '--description-tag', help='Tag that appears in generated tables and plots.', default='DESCRIPTION_TAG')
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('--fdr-threshold', help='FDR threshold for defining directed interactions.', default=0.05)
parser.add_argument('--nominal-alpha-max', help='Maximum nominal alpha for which an FDR is estimated.',
                    default=0.025)
parser.add_argument('--nominal-alpha-step', help='P-value threshold step size used for FDR procedure.',
                    default=0.00001)
parser.add_argument('-n','--iter-num', help='Number of iterations for randomization.', default=100)
parser.add_argument('--random-seed', help='Random seed used for FDR procedure.',
                    default=None)
parser.add_argument('-t','--thread-num', help='Number of threads used for randomization.', default=0)
parser.add_argument('--p-value-threshold', help='P-value threshold for defining directed interactions.', default=None)
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.',
                    required=False)

args = parser.parse_args()
out_prefix = args.out_prefix
description_tag = args.description_tag
diachromatic_interaction_file = args.diachromatic_interaction_file
fdr_threshold = float(args.fdr_threshold)
nominal_alpha_max = float(args.nominal_alpha_max)
nominal_alpha_step = float(args.nominal_alpha_step)
iter_num = int(args.iter_num)
if args.random_seed is None:
    random_seed = args.random_seed
else:
    random_seed = int(args.random_seed)
thread_num = int(args.thread_num)
p_value_threshold = args.p_value_threshold
enriched_digests_file = args.enriched_digests_file

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --description-tag: " + description_tag + '\n'
parameter_info += "\t[INFO] --diachromatic-interaction-file:" + '\n'
parameter_info += "\t\t[INFO] " + diachromatic_interaction_file + '\n'
parameter_info += "\t[INFO] --p-value-threshold: " + str(p_value_threshold) + '\n'
if args.p_value_threshold is None:
    parameter_info += "\t\t[INFO] Will determine a P-value threshold so that the FDR is kept below: " + str(
        fdr_threshold) + '\n'
    parameter_info += "\t\t[INFO] --fdr-threshold: " + "{:.5f}".format(fdr_threshold) + '\n'
    parameter_info += "\t\t[INFO] --nominal-alpha-max: " + "{:.5f}".format(nominal_alpha_max) + '\n'
    parameter_info += "\t\t[INFO] --nominal-alpha-step: " + "{:.5f}".format(nominal_alpha_step) + '\n'
    parameter_info += "\t\t[INFO] --iter-num: " + str(iter_num) + '\n'
    parameter_info += "\t\t[INFO] --random-seed: " + str(random_seed) + '\n'
    parameter_info += "\t\t[INFO] --thread-num: " + str(thread_num) + '\n'
    parameter_info += "\t\t[INFO] Use '--fdr-threshold' to set your own FDR threshold." + '\n'
    parameter_info += "\t\t[INFO] Or use '--p-value-threshold' to skip the FDR procedure." + '\n'
else:
    p_value_threshold = float(p_value_threshold)
    parameter_info += "\t\t[INFO] Will use this P-value threshold instead of the one determined by the FDR procedure." + '\n'
    parameter_info += "\t\t[INFO] We use the negative of the natural logarithm of the P-values." + '\n'
    parameter_info += "\t\t\t[INFO] The chosen threshold corresponds to: -ln(" + str(p_value_threshold) + ") = " + str(
        -log(p_value_threshold)) + '\n'

if enriched_digests_file is not None:
    parameter_info += "\t[INFO] --enriched-digests-file: " + enriched_digests_file

print(parameter_info)

### Perform analysis
####################

# Load interactions
interaction_set = DiachromaticInteractionSet()
# To save memory, we only read interactions that can be significant at 0.05.
min_rp_num, min_rp_num_pval = interaction_set._p_values.find_smallest_significant_n(0.05)
interaction_set.parse_file(diachromatic_interaction_file, min_rp_num=min_rp_num, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()
print()

# Determine P-value threshold so that the FDR is kept below a threshold
if p_value_threshold is None:

    # Create list of nominal alphas
    nominal_alphas = arange(nominal_alpha_step, nominal_alpha_max + nominal_alpha_step, nominal_alpha_step)
    if 0.01 not in nominal_alphas:
        nominal_alphas.append(0.01)

    # Perform randomization procedure
    randomize_fdr = RandomizeInteractionSet(random_seed=random_seed)
    fdr_info_dict = randomize_fdr.perform_randomization_analysis(
        interaction_set=interaction_set,
        nominal_alphas=nominal_alphas,
        iter_num=iter_num,
        thread_num=thread_num,
        verbose=True)
    print()

    # Determine P-value threshold
    result_index = randomize_fdr.get_largest_nominal_alpha_index_at_chosen_fdr_thresh(fdr_threshold)
    p_value_threshold = fdr_info_dict['RESULTS']['NOMINAL_ALPHA'][result_index]

    # Get randomization report for the determined P-value threshold
    fdr_info_info_report = randomize_fdr.get_randomization_info_report(p_value_threshold)

    # Get table row for randomization for the determined P-value threshold
    fdr_info_info_table_row = randomize_fdr.get_randomization_info_table_row(
        nominal_alphas_selected = [p_value_threshold, 0.01],
        description=description_tag.replace(' ','_'))

    # Get entire table with randomization results
    fdr_info_info_table = randomize_fdr.get_randomization_info_table_row(
        description=description_tag.replace(' ','_'))

    # Create plot with Z-score and FDR for each nominal alpha
    randomize_fdr.get_randomization_info_plot_at_chosen_fdr_threshold(
        chosen_fdr_threshold = fdr_threshold,
        pdf_file_name = out_prefix + "_randomization_plot_fdr.pdf",
        description = description_tag)

    # Create randomization histogram for the determined P-value threshold
    randomize_fdr.get_randomization_info_plot(
        nominal_alpha_selected = p_value_threshold,
        pdf_file_name  = out_prefix + "_randomization_histogram_at_threshold.pdf",
        description = description_tag + " - At determined P-value threshold")

    # Create randomization histogram for a nominal alpha of 0.01
    randomize_fdr.get_randomization_info_plot(
        nominal_alpha_selected = 0.01,
        pdf_file_name  = out_prefix + "_randomization_histogram_at_001.pdf",
        description = description_tag + " - At a nominal alpha of 0.01")

# Calculate P-values and assign interactions to 'DI' or 'UI'
interaction_set.evaluate_and_categorize_interactions(p_value_threshold, verbose=True)
eval_cat_info_report = interaction_set.get_eval_cat_info_report()
eval_cat_info_table_row = interaction_set.get_eval_cat_info_table_row(out_prefix)
print()

# Select undirected reference interactions from 'UI'
interaction_set.select_reference_interactions(verbose=True)
select_ref_info_report = interaction_set.get_select_ref_info_report()
select_ref_info_table_row = interaction_set.get_select_ref_info_table_row(out_prefix)
print()

# Write Diachromatic interaction file with two additional columns for P-value and interaction category
f_name_interactions = out_prefix + "_evaluated_and_categorized_interactions.tsv.gz"
interaction_set.write_diachromatic_interaction_file(target_file=f_name_interactions, verbose=True)
write_file_info_report = interaction_set.get_write_file_info_report()
print()

### Create file with summary statistics
#######################################


f_name_summary = out_prefix + "_evaluated_and_categorized_summary.txt"
out_fh_summary = open(f_name_summary, 'wt')

# Chosen parameters
out_fh_summary.write(parameter_info + '\n')

# Report on reading files
out_fh_summary.write(read_file_info_report + '\n')

# Report on the determination of the P-value threshold using the FDR procedure
if args.p_value_threshold is None:
    out_fh_summary.write(fdr_info_info_report + '\n')
    out_fh_summary.write(fdr_info_info_table_row + '\n')

    # Write entire randomization table to file
    f_name_fdr_info_info_table = out_prefix + "_randomization_table.txt"
    out_fh_table = open(f_name_fdr_info_info_table, 'wt')
    out_fh_table.write(fdr_info_info_table + '\n')
    out_fh_table.close()

# Report on evaluation and categorization interactions
out_fh_summary.write(eval_cat_info_report + '\n')
out_fh_summary.write(eval_cat_info_table_row + '\n')

# Report on selection of reference interactions
out_fh_summary.write(select_ref_info_report + '\n')
out_fh_summary.write(select_ref_info_table_row + '\n')

# Report on writing the file
out_fh_summary.write(write_file_info_report + '\n')

# Report on generated files
generated_file_info = "[INFO] Generated files:" + '\n'
generated_file_info += "\t[INFO] " + f_name_summary + '\n'
generated_file_info += "\t[INFO] " + f_name_interactions + '\n'
if args.p_value_threshold is None:
    generated_file_info += "\t[INFO] " + f_name_fdr_info_info_table + '\n'
    generated_file_info += "\t[INFO] " + out_prefix + "_randomization_plot_fdr.pdf" + '\n'
    generated_file_info += "\t[INFO] " + out_prefix + "_randomization_plot_threshold.pdf" + '\n'
    generated_file_info += "\t[INFO] " + out_prefix + "_randomization_plot_001.pdf" + '\n'
out_fh_summary.write(generated_file_info)
out_fh_summary.close()

print(generated_file_info + '\n')
