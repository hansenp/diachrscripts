import argparse
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

# Parse command line
####################

parser = argparse.ArgumentParser(description='Re-select reference interactions.')
parser.add_argument('-o', '--out-prefix',
                    help='Common prefix for all generated files, which can also contain the path.',
                    default='OUT_PREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file',
                    help='Input file in Diachromatic interaction format.',
                    required=True)

# Read command line arguments into variables
args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.diachromatic_interaction_file

# Read interaction file to interaction set
interaction_set = DiachromaticInteractionSet(rpc_rule = 'ht')
interaction_set.parse_file(i_file=diachromatic_interaction_file,
                           verbose=True)

# Re-select reference
interaction_set.select_reference_interactions_2x(verbose=True)

# Write interaction file
interaction_set.write_diachromatic_interaction_file(target_file=out_prefix + "_re_selected_interactions.tsv.gz",
                                                    verbose=True)
