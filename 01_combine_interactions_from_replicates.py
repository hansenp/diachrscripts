#!/usr/bin/env python

"""
This script takes a path to a directory containing Diachromatic interaction files and combines interactions that occur
in a specified number of files into one interaction with summed simple and twisted read pair counts.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
import os
from diachr.diachromatic_interaction_parser import DiachromaticInteractionParser


### Parse command line
######################

parser = argparse.ArgumentParser(description='Combine interactions that occur in a specified number of replicates.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i','--interaction-files-path', help='Path to directory with Diachromatic interaction files', required=True)
parser.add_argument('-r','--required-replicates', help='Required number of replicates.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_files_path = args.interaction_files_path
required_replicates = int(args.required_replicates)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --interaction-files-path: " + interaction_files_path)
print("\t[INFO] --required-replicates: " + str(required_replicates))


### Iterate and combine interactions from different replicates
##############################################################

def get_gzip_tsv_files(dir):
    """
    Get list of all gzip files in a directory
    """
    gzfiles = []
    for file in os.listdir(dir):
        if file.endswith(".tsv.gz"):
            gzpath = os.path.join(dir, file)
            gzfiles.append(gzpath)
    return gzfiles


# Get list of interaction files under given path
gz_files = get_gzip_tsv_files(interaction_files_path)
if len(gz_files) < int(required_replicates):
    print("[FATAL] Not enough replicates. Must be at least " + str(required_replicates) + " But there are only " + str(len(gz_files)) + " files.")
    exit(1)

# Get list of Diachromatic interaction objects
interaction_set = DiachromaticInteractionParser()
print("[INFO] Will parse all gz files in " + interaction_files_path)
for gz_file in gz_files:
    interaction_set.parse_file(gz_file)

# Print information about read files
print(interaction_set.get_file_dict_info())

# Write interactions that occur in the required number of replicates to file
f_name = out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_interactions.tsv.gz"
write_info = interaction_set.write_diachromatic_interaction_file(target_file_name=f_name, required_replicates=required_replicates)

# Print information about written interactions
print(write_info)

# Create file that contains all the information
f_name =  out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_summary.txt"
out_fh = open(f_name, 'wt')
out_fh.write(interaction_set.get_file_dict_info() + '\n')
out_fh.write(write_info + '\n')
out_fh.close()

print(len(interaction_set.interaction_list))


n_in
for i in interaction_set.interaction_list:
    print(i.has_data_for_required_replicate_num(2))
    print(i.chrA + "\t" + str(i.fromA) + "\t" + str(i.toA) + "\t" + i.chrA + "\t")
