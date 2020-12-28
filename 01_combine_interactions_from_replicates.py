#!/usr/bin/env python

"""
This script takes a path to a directory containing Diachromatic interaction files and combines interactions that occur
in a specified number of files into one interaction with summed simple and twisted read pair counts.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
import os
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet


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

parameter_info = "[INFO] " + "Input parameters" + '\n'
parameter_info += "\t[INFO] --out-prefix: " + out_prefix + '\n'
parameter_info += "\t[INFO] --interaction-files-path: " + interaction_files_path + '\n'
parameter_info += "\t[INFO] --required-replicates: " + str(required_replicates) + '\n'
print(parameter_info)

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

# Read interaction files
interaction_set = DiachromaticInteractionSet()
read_info = "[INFO] Will parse all gz files in " + interaction_files_path + '\n'
for gz_file in gz_files:
    interaction_set.parse_file(gz_file)

# Print information about read files
file_info_dict = interaction_set.get_file_info_dict()
file_num = len(file_info_dict['I_FILE']) - 1
read_info += "\t[INFO] Read interaction data from " + str(file_num) + " files:" + '\n'
for i in range(0, file_num):
    read_info += "\t\t[INFO] " + str(file_info_dict['I_NUM'][i]) + " interactions from " + file_info_dict['I_FILE'][i] + '\n'
read_info += "\t[INFO] The union of all interactions has " + str(file_info_dict['I_NUM'][file_num]) + " interactions." + '\n'
read_info += "[INFO] ... done." + '\n'
print(read_info)

# Write interactions that occur in the required number of replicates to file
f_name = out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_interactions.tsv.gz"
write_info_dict = interaction_set.write_diachromatic_interaction_file(target_file=f_name, required_replicates=required_replicates)

# Print information about written interactions
write_info = "[INFO] Writing interactions that occur in at least " + str(required_replicates) + " replicates to: " + str(write_info_dict['TARGET_FILE'][0]) + '\n'
write_info += "\t[INFO] Interactions that occur in at least " + str(write_info_dict['REQUIRED_REPLICATES'][0]) + " replicates: " + str(write_info_dict['HAS_ALL_DATA'][0]) + '\n'
write_info += "\t[INFO] Other interactions: " + str(write_info_dict['INCOMPLETE_DATA'][0]) + '\n'
write_info += "[INFO] ... done." + '\n'
print(write_info)

# Create a table row with all information
table_row = "TARGET_FILE" + "\t" + "I_NUMS" + "\t" + "REQUIRED_REPLICATES" + "\t" + "HAS_ALL_DATA" + "\t" + "INCOMPLETE_DATA" + '\n'
table_row += str(write_info_dict['TARGET_FILE'][0]) + "\t" +\
              str(file_info_dict['I_NUM'][:-1]) + "\t" +\
              str(write_info_dict['REQUIRED_REPLICATES'][0]) + "\t" +\
              str(write_info_dict['HAS_ALL_DATA'][0]) + "\t" +\
              str(write_info_dict['INCOMPLETE_DATA'][0]) + '\n'
print(table_row)

# Create file that contains all the information
f_name =  out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_summary.txt"
out_fh = open(f_name, 'wt')
out_fh.write(parameter_info + '\n')
out_fh.write(read_info + '\n')
out_fh.write(write_info + '\n')
out_fh.write(table_row)
out_fh.close()



