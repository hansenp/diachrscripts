"""
This script takes a file in Diachromatic interaction format and a P-value threshold, calculates the P-values of
interactions, assigns interactions to directed (DI) and undirected interactions (UI), selects undirected reference
interactions (UIR) from UI and creates a new file with two additional columns for P-value and interaction category.

You can find documentation on this script in the relevant section in the RTD of this repository.

The individual steps that are carried out in this script are demonstrated in the following Jupyter notebook:

       diachrscripts/jupyter_notebooks/evaluate_and_categorize_interactions.ipynb
"""

import argparse
from numpy import log
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


### Parse command line
######################

parser = argparse.ArgumentParser(description='Evaluate and categorize interactions and select unndirected reference interactions.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('--p-value-threshold', help='P-value threshold for directed interactions.', default=0.01)
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.', required=False)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.diachromatic_interaction_file
p_value_threshold = float(args.p_value_threshold)
enriched_digests_file = args.enriched_digests_file

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file + '\n'
parameter_info += "\t[INFO] --p-value-threshold: " + str(p_value_threshold) + '\n'
if enriched_digests_file != None:
    parameter_info += "\t[INFO] --enriched-digests-file: " + enriched_digests_file

parameter_info += "\t[INFO] We use the negative of the natural logarithm of the P-values." + '\n'
parameter_info += "\t\t[INFO] The chosen threshold corresponds to: -ln(" + str(p_value_threshold) + ") = " + str(-log(p_value_threshold)) + '\n'

print(parameter_info)


### Perform analysis
####################

# Load interactions
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()

# Calculate P-values and assign interactions to 'DI' or 'UI'
interaction_set.evaluate_and_categorize_interactions(p_value_threshold, verbose=True)
eval_cat_info_report = interaction_set.get_eval_cat_info_report()
eval_cat_info_table_row = interaction_set.get_eval_cat_info_table_row(out_prefix)

# Select undirected reference interactions from 'UI'
interaction_set.select_reference_interactions(verbose=True)
select_ref_info_report = interaction_set.get_select_ref_info_report()
select_ref_info_table_row = interaction_set.get_select_ref_info_table_row(out_prefix)

# Write Diachromatic interaction file with two additional columns
f_name_interactions = out_prefix + "_evaluated_and_categorized_interactions.tsv.gz"
interaction_set.write_diachromatic_interaction_file(target_file=f_name_interactions, verbose=True)
write_file_info_report = interaction_set.get_write_file_info_report()


### Create file with summary statistics
#######################################

f_name_summary = out_prefix + "_evaluated_and_categorized_summary.txt"
out_fh = open(f_name_summary, 'wt')

# Chosen parameters
out_fh.write(parameter_info + '\n')

# Report on reading files
out_fh.write(read_file_info_report + '\n')

# Report on evaluation and categorization interactions
out_fh.write(eval_cat_info_report + '\n')
out_fh.write(eval_cat_info_table_row + '\n')

# Report on selection of reference interactions
out_fh.write(select_ref_info_report + '\n')
out_fh.write(select_ref_info_table_row + '\n')

# Report on writing the file
out_fh.write(write_file_info_report + '\n')

# Report on generated files
generated_file_info = "[INFO] Generated files:" + '\n'
generated_file_info += "\t[INFO] " + f_name_summary + '\n'
generated_file_info += "\t[INFO] " + f_name_interactions + '\n'
out_fh.write(generated_file_info)

out_fh.close()

print(generated_file_info)
