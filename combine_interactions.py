#!/usr/bin/env python
import argparse
import os
import gzip
from collections import defaultdict


### Parse command line
######################

parser = argparse.ArgumentParser(description='Discard all interactions that are not significant in a given number of replicates.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--out-dir', help='Directory for output.', default='.')
parser.add_argument('--interaction-files-path', help='Path to directory with gt1 gzip files')
parser.add_argument('--n-experiments', help='Path to directory with gt1 gzip files')

args = parser.parse_args()
out_prefix = args.out_prefix
out_dir = args.out_dir
interaction_files_path = args.interaction_files_path
n_experiments = int(args.n_experiments)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Path to gzipped interaction files: " + interaction_files_path)


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
        self.chrA = F[0]
        self.fromA = F[1]
        self.toA = F[2]
        self.statusA = F[3]
        self.chrB = F[4]
        self.fromB = F[5]
        self.toB = F[6]
        self.statusB = F[7]
        self.twisted = []
        self.simple = []
        self.counts_string = []
        counts = F[8]  # a string such as "1:3"
        self._append_twisted_and_simple_counts(counts)

    def append_interaction_data(self, F):
        """
        This method is called if we already have one or more observations
        for this interaction. We have thus compared the key already.
        We can just add the interaction counts.
        """
        counts = F[8]  # a string such as "1:3"
        self._append_twisted_and_simple_counts(counts)

    def has_data_for_all_experiments(self, k):
        return len(self.twisted) == k

    def _append_twisted_and_simple_counts(self, counts):
        """
        'private' method to add the counts
        """
        cc = counts.split(":")

        if len(cc) != 2:
            raise TypeError("Malformed counts string {}".format(cc))
        self.simple.append(int(cc[0]))
        self.twisted.append(int(cc[1]))
        self.counts_string.append(counts)

    @staticmethod
    def make_key(F):
        return "{}:{}-{}_{}:{}-{}".format(F[0], F[1], F[2], F[4], F[5], F[6])

    def output_summary_2(self):
        cstring = ",".join(self.counts_string)
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}\t{}".format(self.chrA,
                                                                              self.fromA,
                                                                              self.toA,
                                                                              self.statusA,
                                                                              self.chrB,
                                                                              self.fromB,
                                                                              self.toB,
                                                                              self.statusB,
                                                                              sum(self.simple),
                                                                              sum(self.twisted),
                                                                              cstring)

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
            if len(F) != 9:
                raise TypeError("Malformed line {}".format(line))
            mykey = Interaction.make_key(F)
            if mykey in interaction_dict:
                iaction = interaction_dict[mykey]
                iaction.append_interaction_data(F)
            else:
                interaction_dict[mykey] = Interaction(F)

    return n_iteraction


### Iterate interactions
########################

print("[INFO] Will parse all gz files in", interaction_files_path)

d = defaultdict(Interaction)

n_interactions = []
gzfiles = get_gzip_tsv_files(interaction_files_path)
for f in gzfiles:
    print("\t[INFO] Extracting interactions from " + f + ".")
    n = parse_gzip_tsv_file(f, d)
    n_interactions.append(n)
if len(gzfiles) == 0:
    print("[FATAL] Did not find any gzipped files. Note: Need to give path to directory with *.tsv.gz files.")
    exit(1)

print("\t[INFO] The union of all interactions has " + str(len(d)) + " interactions.")

if len(gzfiles) < int(n_experiments):
    print("[FATAL] Not enough replicates. Must be at least " + str(n_experiments) + " But there are only " + str(len(gzfiles)) + " files.")
    exit(1)


fname = out_dir + "/" + out_prefix + "_alt_interactions.txt.gz"
outfh = gzip.open(fname, 'wt')
print("[INFO] Writing interactions to file ...")
n_has_all_data = 0
n_incomplete_data = 0
for key, iaction in d.items():
    if not iaction.has_data_for_all_experiments(n_experiments):
        n_incomplete_data += 1
    else:
        n_has_all_data += 1
        outfh.write(iaction.output_summary_2() + "\n")

    if (n_incomplete_data + n_has_all_data)%1000000==0:
        print("\t[INFO] " + (str(n_incomplete_data + n_has_all_data)) + " interactions processed ...")
outfh.close()

fname = out_dir + "/" + out_prefix + "_alt_summary.txt"
outfh = open(fname, 'wt')
outfh.write(str(out_prefix) + "\t" + str(n_interactions) + "\t" + str(n_experiments) + "\t" + str(n_has_all_data) + "\t" + str(n_incomplete_data) + "\n")
outfh.close()

print("Total number of interactions: " + str(n_interactions))
print("[INFO] I nteractions with {} data points: {}, lacking interactions: {}".format(n_experiments, n_has_all_data, n_incomplete_data))
print("[INFO] We wrote all interactions to file: {}".format(fname))
print("[INFO] ... done.")
