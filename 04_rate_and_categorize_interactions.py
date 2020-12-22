"""
This script takes a file in Diachromatic interaction format and a P-value threshold, calculates the P-values of
interactions, assigns interactions to directed (DI) and undirected interactions (UI), selects undirected reference
interactions (UIR) from UI and creates a new file with two additional columns for P-value and interaction category.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
from numpy import log
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


### Parse command line
######################

parser = argparse.ArgumentParser(description='Rate and categorize interactions and select unndirected reference interactions.')
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

# Transform P-value threshold
nln_p_value_threshold = -log(p_value_threshold)

# Load interactions
interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file)

# Calculate P-values and assign interactions to 'DI' or 'UI'
rate_and_cat_report, rate_and_cat_summary_stat_dict = interaction_set.rate_and_categorize_interactions(nln_p_value_threshold, verbose=True)

# Select undirected reference interactions from 'UI'
select_ref_report, select_ref_summary_stat_dict, missing_ref_info = interaction_set.select_reference_interactions(verbose=True)

# Write Diachromatic interaction file with two additional columns
f_name = out_prefix + "_rated_and_categorized_interactions.tsv.gz"
interaction_set.write_diachromatic_interaction_file(target_file_name=f_name, verbose=True)


### Create file with summary statistics
#######################################

f_name = out_prefix + "_rate_and_categorize_summary.txt"
out_fh = open(f_name, 'wt')

# Chosen parameters
out_fh.write(parameter_info + '\n')

out_fh.write(rate_and_cat_report + '\n\n')

out_fh.write(
    "OUT_PREFIX" + '\t' +
    "NLN_PVAL_THRESH" + '\t' +
    "MIN_RP" + '\t' +
    "MIN_RP_NLN_PVAL" + '\t' +
    "N_PROCESSED" + '\t' +
    "N_DISCARDED" + '\t' +
    "N_UNDIRECTED" + '\t'
    "N_DIRECTED" + '\n'
)
out_fh.write(
    out_prefix + '\t' +
    "{:.2f}".format(rate_and_cat_summary_stat_dict['NLN_PVAL_THRESH']) + '\t' +
    str(rate_and_cat_summary_stat_dict['MIN_RP']) + '\t' +
    "{:.2f}".format(-log(rate_and_cat_summary_stat_dict['MIN_RP_PVAL'])) + '\t' +
    str(rate_and_cat_summary_stat_dict['N_PROCESSED']) + '\t' +
    str(rate_and_cat_summary_stat_dict['N_DISCARDED']) + '\t' +
    str(rate_and_cat_summary_stat_dict['N_UNDIRECTED']) + '\t' +
    str(rate_and_cat_summary_stat_dict['N_DIRECTED']) + '\n\n'
)
out_fh.write(select_ref_report + '\n\n')

out_fh.write(
    "OUT_PREFIX" + '\t' +

    "DI_NN" + '\t' +
    "DI_NE" + '\t' +
    "DI_EN" + '\t' +
    "DI_EE" + '\t' +

    "UIR_NN" + '\t' +
    "UIR_NE" + '\t' +
    "UIR_EN" + '\t' +
    "UIR_EE" + '\t' +

    "M_UIR_NN" + '\t' +
    "M_UIR_NE" + '\t' +
    "M_UIR_EN" + '\t' +
    "M_UIR_EE" + '\t' +

    "UI_NN" + '\t' +
    "UI_NE" + '\t' +
    "UI_EN" + '\t' +
    "UI_EE" + '\n'
)

out_fh.write(
    out_prefix + '\t' +

    str(select_ref_summary_stat_dict["DI_NN"]) + '\t' +
    str(select_ref_summary_stat_dict["DI_NE"]) + '\t' +
    str(select_ref_summary_stat_dict["DI_EN"]) + '\t' +
    str(select_ref_summary_stat_dict["DI_EE"]) + '\t' +

    str(select_ref_summary_stat_dict["UIR_NN"]) + '\t' +
    str(select_ref_summary_stat_dict["UIR_NE"]) + '\t' +
    str(select_ref_summary_stat_dict["UIR_EN"]) + '\t' +
    str(select_ref_summary_stat_dict["UIR_EE"]) + '\t' +

    str(select_ref_summary_stat_dict["M_UIR_NN"]) + '\t' +
    str(select_ref_summary_stat_dict["M_UIR_NE"]) + '\t' +
    str(select_ref_summary_stat_dict["M_UIR_EN"]) + '\t' +
    str(select_ref_summary_stat_dict["M_UIR_EE"]) + '\t' +

    str(select_ref_summary_stat_dict["UI_NN"]) + '\t' +
    str(select_ref_summary_stat_dict["UI_NE"]) + '\t' +
    str(select_ref_summary_stat_dict["UI_EN"]) + '\t' +
    str(select_ref_summary_stat_dict["UI_EE"]) + '\n'
)

out_fh.close()
