"""
In this script, the P-values for directionality of interactions are calculated and, based on these P-values,
the interactions are categorized into directed (DI) und undirected interactions (UI).
Undirected reference interactions (UIR) are then selected from UI, which are comparable to DI with respect to
to the total number of read pairs (n).

The input consists of a file in Diachromatic interaction format and a P-value threshold that was determined using our
FDR procedure. The output is again a file in Diachromatic interaction format, but there are two additional columns
on the right for the calculated P-values and interaction categories.
"""

import argparse
from numpy import log
from diachr.diachromatic_interaction_parser import DiachromaticInteractionParser

### Parse command line
######################

parser = argparse.ArgumentParser(description='Rate and categorize interactions and select unndirected reference interactions.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file', help='Diachromatic interaction file.', required=True)
parser.add_argument('--p-value-threshold', help='P-value threshold for directed interactions.', default=0.001)
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.diachromatic_interaction_file
p_value_threshold = float(args.p_value_threshold)
enriched_digests_file = args.enriched_digests_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --p-value-threshold: " + str(p_value_threshold))
print("\t[INFO] --enriched-digests-file: " + enriched_digests_file)

# Transform P-value threshold
nln_p_value_threshold = -log(p_value_threshold)
print("[INFO] We use the negative of the natural logarithm of the P-values.")
print("\t[INFO] The chosen threshold corresponds to: -ln(" + str(p_value_threshold) + ") = " + str(nln_p_value_threshold))

# Load interactions with 'DiachromaticInteractionParser'
interaction_set = DiachromaticInteractionParser(enriched_digests_file)
interaction_set.parse_file(diachromatic_interaction_file)

# Calculate P-values and assign interactions to 'DI' or 'UIR'
interaction_set.rate_and_categorize_interactions(nln_p_value_threshold)


interaction_set.select_reference_interactions()

f_name = out_prefix + "_rated_and_categorized_interactions" + ".tsv.gz"
interaction_set.write_diachromatic_interaction_file(target_file_name=f_name)




# FIRST PASS: CALCULATE P-VALUES AND DEFINE DIRECTED AND UNDIRECTED INTERACTIONS

# SECOND PASS: SELECT UNDIRECTED REFERENCE INTERACTIONS

# WRITE FILE IN DIACHROMATIC INTERACTION FORMAT WITH TWO ADDITIONAL COLUMNS


