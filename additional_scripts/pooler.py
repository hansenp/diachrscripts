#!/usr/bin/env python

"""
This script takes a path to a directory containing Diachromatic interaction files and combines interactions that occur
in a specified number of files into one interaction with summed simple and twisted read pair counts.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
import os
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

# Parse command line
####################

parser = argparse.ArgumentParser(description='Combine interactions that occur in a specified number of replicates.')
parser.add_argument('-o', '--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i', '--interaction-files-path', help='Path to directory with Diachromatic interaction files',
                    required=True)
parser.add_argument('-r', '--required-replicates', help='Required number of replicates.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_files_path = args.interaction_files_path
required_replicates = int(args.required_replicates)

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --interaction-files-path: " + interaction_files_path + '\n'
parameter_info += "\t[INFO] --required-replicates: " + str(required_replicates) + '\n'

print(parameter_info)


# Get list of interaction files under given path
################################################

def get_gzip_tsv_files(path):
    """
    Get list of all gzip files in a directory
    """
    gz_files = []
    for file in os.listdir(path):
        if file.endswith(".tsv.gz"):
            gz_path = os.path.join(path, file)
            gz_files.append(gz_path)
    return gz_files


gz_files = get_gzip_tsv_files(interaction_files_path)
if len(gz_files) < int(required_replicates):
    print("[FATAL] Not enough replicates. Must be at least " + str(required_replicates) + " But there are only " + str(
        len(gz_files)) + " files.")
    exit(1)

# Perform analysis
##################

# Read interaction files
interaction_set = DiachromaticInteractionSet()
for gz_file in gz_files:
    interaction_set.parse_file(i_file=gz_file, verbose=True)
read_file_info_report = interaction_set.get_read_file_info_report()
read_file_info_table_row = interaction_set.get_read_file_info_table_row()
print()

# Write interactions that occur in the required number of replicates to file
f_name_interactions = out_prefix + "_at_least_" + str(required_replicates) + "_combined_interactions.tsv.gz"
interaction_set.write_diachromatic_interaction_file(target_file=f_name_interactions,
                                                    required_replicates=required_replicates, verbose=True)
write_file_info_report = interaction_set.get_write_file_info_report()
write_file_info_table_row = interaction_set.get_write_file_info_table_row()

# Create file with summary statistics
#####################################

f_name_summary = out_prefix + "_at_least_" + str(required_replicates) + "_combined_summary.txt"

out_fh = open(f_name_summary, 'wt')

# Chosen parameters
out_fh.write(parameter_info + '\n')

# Report on reading files
out_fh.write(read_file_info_report + '\n')
out_fh.write(read_file_info_table_row + '\n')

# Report on writing the file
out_fh.write(write_file_info_report + '\n')
out_fh.write(write_file_info_table_row + '\n')

# Report on generated files
generated_file_info = "[INFO] Generated files:" + '\n'
generated_file_info += "\t[INFO] " + f_name_summary + '\n'
generated_file_info += "\t[INFO] " + f_name_interactions + '\n'
out_fh.write(generated_file_info)

out_fh.close()
