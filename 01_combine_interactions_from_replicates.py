#!/usr/bin/env python

"""
This takes a path to a directory containing Diachromatic interaction files and combines interactions that occur in a
specified number of files into one interaction with summed simple and twisted read pair counts.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse
from diachr.diachromatic_parser import DiachromaticParser
import os
import gzip
from collections import defaultdict
from diachr.diachromatic_interaction import DiachromaticInteraction
from diachr.diachromatic_interaction_parser import DiachromaticInteractionParser


### Parse command line
######################

parser = argparse.ArgumentParser(description='Combine interactions that occur in a specified number of replicates.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('-i','--interaction-files-path', help='Path to directory with gzip files', required=True)
parser.add_argument('-r','--required-replicates', help='Required number of replicates.', required=True)
parser.add_argument('-m','--usemod', help="Use new module", dest='usemod', action='store_true')

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_files_path = args.interaction_files_path
required_replicates = int(args.required_replicates)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --interaction-files-path: " + interaction_files_path)
print("\t[INFO] --required-replicates: " + str(required_replicates))
print("\t[INFO] --usemod: " + str(args.usemod))


############   Uses new module/class ###############
if args.usemod:
    parser = DiachromaticParser(path=interaction_files_path, required_replicates=required_replicates, prefix=out_prefix)
    parser.write_combined_interaction_file()
    exit(1)
############   End ###############


### Define auxiliary classes and methods
########################################

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


class Interaction:
    """
    Class to represent an interaction between
    two different parts of the genome with
    status and count.
    """

    def __init__(self, F):
        """
        @param F an array with 9 elements that represent an interaction
        """
        global chr_id
        global chr_dict

        # map status flag pair to int
        if F[3] == 'I' and F[7] == 'I':
            self.s_flag_int = 0
        elif F[3] == 'I' and F[7] == 'A':
            self.s_flag_int = 1
        elif F[3] == 'A' and F[7] == 'I':
            self.s_flag_int = 2
        else: # AA
            self.s_flag_int = 3

        if F[0] in chr_dict:
            self.chrA = chr_dict[F[0]]
        else:
            self.chrA = chr_id
            chr_dict[F[0]] = chr_id
            chr_id += 1
        self.fromA = F[1]
        self.toA = F[2]
        if F[4] in chr_dict:
            self.chrB = chr_dict[F[4]]
        else:
            self.chrB = chr_id
            chr_dict[F[4]] = chr_id
            chr_id += 1
        self.fromB = F[5]
        self.toB = F[6]
        self.twisted = 0
        self.simple = 0
        self.rep_num = 1
        self._append_twisted_and_simple_counts(F[8]) # F[8] contains a string such as "1:3"

    def append_interaction_data(self, F):
        """
        This method is called if we already have one or more observations
        for this interaction. We have thus compared the key already.
        We can just add the interaction counts.
        """
        self.rep_num += 1
        self._append_twisted_and_simple_counts(F[8]) # F[8] contains a string such as "1:3"

    def has_data_for_required_replicate_num(self, k):
        return k <= self.rep_num

    def _append_twisted_and_simple_counts(self, counts):
        """
        'private' method to add the counts
        """
        cc = counts.split(":")
        if len(cc) != 2:
            raise TypeError("Malformed counts string {}".format(cc))
        self.simple = self.simple + int(cc[0])
        self.twisted = self.twisted + int(cc[1])

    def output_summary_2(self):
        global inv_chr_dict
        if self.s_flag_int == 0:
            statusA = 'I'
            statusB = 'I'
        elif self.s_flag_int == 1:
            statusA = 'I'
            statusB = 'A'
        elif self.s_flag_int == 2:
            statusA = 'A'
            statusB = 'I'
        else:
            statusA = 'A'
            statusB = 'A'
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}".format(inv_chr_dict[self.chrA],
                                                                              self.fromA,
                                                                              self.toA,
                                                                              statusA,
                                                                              inv_chr_dict[self.chrB],
                                                                              self.fromB,
                                                                              self.toB,
                                                                              statusB,
                                                                              self.simple,
                                                                              self.twisted)

    def __hash__(self):
        return hash((self.chrA, self.fromA, self.chrB, self.fromB))

    def __eq__(self, other):
        return self.fromA == other.fromA and self.fromB == other.fromB and self.chrA == other.chrA and self.chrB == other.chrB

def parse_gzip_tsv_file(file, interaction_dict):
    """
    Parse output of diachromatic interaction file
    with goal of getting the interaction counts for
    simple/twisted
    chr5\t169261149\t169264959\tI\tchr5\t169981022\t169985681\tI\t0:2
    """
    n_iteraction = 0
    with gzip.open(file, 'r' + 't') as f:
        for line in f:
            n_iteraction += 1
            F = line.rstrip().split('\t')
            if len(F) < 9:
                raise TypeError("Malformed line {}".format(line))
            iaction = Interaction(F)
            if iaction in interaction_dict:
                 interaction_dict[iaction].append_interaction_data(F)
            else:
                interaction_dict[iaction] = iaction
    f.close()

    return n_iteraction


### Iterate interactions
########################
# THE FOLLOWING COMMENTED LINES DO NOT USE THE FUNCTIONS FROM THIS SCRIPT,
# BUT ONLY THE CLASSES DiachromaticInteraction and DiachromaticInteractionParser.
# WHICH SIMPLIFIES THE CODE.
# HOWEVER, THIS TASK REQUIRES A LOT OF MEMORY ANND IT STILL HAS TO BE TESTED WHETHER THE
# MANY INTERACTION OBJECTS ARE NOT TOO LARGE.
# IF THIS IS NOT THE CASE, WE COULD DELETE EVERYTHING ELSE FROM THIS SCRIPT.

# Get list of interaction files under given path
gz_files = get_gzip_tsv_files(interaction_files_path)
if len(gz_files) < int(required_replicates):
    print("[FATAL] Not enough replicates. Must be at least " + str(required_replicates) + " But there are only " + str(len(gz_files)) + " files.")
    exit(1)

# Get list of Diachromatic interaction objects
interaction_set = DiachromaticInteractionParser()
print("[INFO] Will parse all gz files in " + interaction_files_path)
for gz_file in gz_files:
    #print("\t[INFO] Reading ", gz_file)
    interaction_set.parse_file(gz_file)

# Print information about read files
print(interaction_set.get_file_dict_info())

# Write interactions that occur in enough replicates to file
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

exit(1)

### Iterate interactions
########################
# THIS IS THE OLD IMPLEMENTATION THAT CAN BE DELETED IF THE NEW IMPLEMENTATION DOES NOT
# USE TOO MUCH MEMORY.

print("[INFO] Will parse all gz files in", interaction_files_path)

chr_dict = dict() # global dictionary maps chromosome string to integer, e.g. chr1 to 4
chr_id = 0

d = defaultdict(Interaction)
n_interactions = []
gzfiles = get_gzip_tsv_files(interaction_files_path)
for f in gzfiles:
    print("\t[INFO] Extracting interactions from " + f + " ...")
    n = parse_gzip_tsv_file(f, d)
    n_interactions.append(n)
if len(gzfiles) == 0:
    print("[FATAL] Did not find any gzipped files. Note: Need to give path to directory with *.tsv.gz files.")
    exit(1)

print("\t[INFO] The union of all interactions has " + str(len(d)) + " interactions.")
print("[INFO] Total numbers of interactions: " + str(n_interactions))

if len(gzfiles) < int(required_replicates):
    print("[FATAL] Not enough replicates. Must be at least " + str(required_replicates) + " But there are only " + str(len(gzfiles)) + " files.")
    exit(1)

fname = out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_interactions.tsv.gz"
outfh = gzip.open(fname, 'wt')
print("[INFO] Writing interactions to file ...")
# reverse dictionary with chromosome IDs
inv_chr_dict = {v: k for k, v in chr_dict.items()}
n_has_required_data = 0
n_incomplete_data = 0
for key, iaction in d.items():
    if not iaction.has_data_for_required_replicate_num(required_replicates):
        n_incomplete_data += 1
    else:
        n_has_required_data += 1
        outfh.write(iaction.output_summary_2() + "\n")

    if (n_incomplete_data + n_has_required_data)%1000000==0:
        print("\t[INFO] " + (str(n_incomplete_data + n_has_required_data)) + " interactions processed ...")
outfh.close()
print("[INFO] We wrote all interactions to file: {}".format(fname))

fname = out_prefix + "_at_least_in_" + str(required_replicates) + "_replicates_summary.txt"
outfh = open(fname, 'wt')
outfh.write("OUT_PREFIX" + "\t" + "INTERACTIONS_NUMBERS" + "\t" + "REQUIRED_INTERACTIONS" + "\t" + "HAS_ALL_DATA" + "\t" + "INCOMPLETE_DATA" + "\n")
outfh.write(str(out_prefix) + "\t" + str(n_interactions) + "\t" + str(required_replicates) + "\t" + str(n_has_required_data) + "\t" + str(n_incomplete_data) + "\n")
outfh.close()

print("[INFO] Interactions with at least {} data points: {}, lacking interactions: {}".format(required_replicates, n_has_required_data, n_incomplete_data))
print("[INFO] We wrote summary statistics to file: {}".format(fname))
print("[INFO] ... done.")
