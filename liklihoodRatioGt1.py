#!/usr/bin/env python
import argparse
import os
import gzip
import numpy as np
from collections import defaultdict
from scipy.stats.distributions import chi2
from scipy.stats import binom
from math import log

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

    def output_summary(self, LR, p, threshold):
        cstring = ",".join(self.counts_string)
        twistedsum = sum(self.twisted)
        simplesum = sum(self.simple)
        if p >= threshold:
            status = "U"
        elif twistedsum > simplesum:
            status = "T"
        elif twistedsum < simplesum:
            status = "S"
        else:
            raise TypeError("Should never happen -- significant p value but T == S!")
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chrA, self.fromA, self.toA,
                                                                               self.statusA,
                                                                               self.chrB, self.fromB, self.toB,
                                                                               self.statusB,
                                                                               sum(self.twisted), sum(self.simple),
                                                                               cstring, LR, p, status)


def parse_gzip_tsv_file(file, interaction_dict):
    """
    Parse output of diachromatic interaction file
    with goal of getting the interaction counts for
    simple/twisted
    chr5\t169261149\t169264959\tI\tchr5\t169981022\t169985681\tI\t0:2
    """
    with gzip.open(file, 'r' + 't') as f:
        for line in f:
            F = line.rstrip().split('\t')
            if len(F) != 9:
                raise TypeError("Malformed line {}".format(line))
            ## If we have already seen this interaction, append the
            ## counts, otherwise create a new Interaction object
            mykey = Interaction.make_key(F)
            if mykey in interaction_dict:
                iaction = interaction_dict[mykey]
                iaction.append_interaction_data(F)
            else:
                interaction_dict[mykey] = Interaction(F)


def likelihood_ratio_test(twisted, simple):
    p1 = 0.5
    tw = sum(twisted)
    si = sum(simple)
    # p2 is estimated by maximum likelihood here:
    p2 = float(tw) / float(tw + si)
    b1 = binom.pmf(tw, (tw + si), p1)
    b2 = binom.pmf(tw, (tw + si), p2)
    if b1 <= 0 or b2 <= 0:
        print b1
        print b2
        print p1
        print p2
        print tw
        print si
        return 100000.00, 2.0
    LR = 2 * (log(b2) - log(b1))
    p = chi2.sf(LR, 1)  # one degree of freedom, chi2 distributed
    return LR, p


### Start Execution

parser = argparse.ArgumentParser(description='Process gt1 interaction files.')
parser.add_argument('--gzdir', help='path to directory with gt1 gzip files')
parser.add_argument('--outdir', help='path to directory to which the results will be written')
parser.add_argument('--outprefix', help='prefix for output files')
args = parser.parse_args()
gzpath = args.gzdir
outdir = args.outdir
outprefix = args.outprefix

print("[INFO] Will parse all gz files in", gzpath)

d = defaultdict(Interaction)

gzfiles = get_gzip_tsv_files(gzpath)
for f in gzfiles:
    print("[INFO] Extracting interactions from ", f)
    parse_gzip_tsv_file(f, d)
if len(gzfiles) == 0:
    print("[FATAL] Did not find any gzipped files. Note: Need to give path to directory with *.tsv.gz files.")
    exit(1)

## When we get here, d has interaction objects for all observed interactions
## Some of them have observations for all experiments

n_experiments = len(gzfiles)
n_has_all_data = 0
n_incomplete_data = 0

i = 0

fname = outdir + "/" + outprefix + "_interactions.txt"
outfh = open(fname, 'w')

for key, iaction in d.items():
    if not iaction.has_data_for_all_experiments(n_experiments):
        n_incomplete_data += 1
    else:
        n_has_all_data += 1
        LR, p = likelihood_ratio_test(iaction.twisted, iaction.simple)
        i += 1
        pthreshold = 0.0001
        outfh.write(iaction.output_summary(LR, p, pthreshold) + "\n")
        #if p < 0.0001:
            #print("{}".format(iaction.output_summary(LR, p, pthreshold)))

print("[INFO] Interactions with all {} data points: {}, lacking data {}".format(n_experiments, n_has_all_data,
                                                                                n_incomplete_data))
print("[INFO] We wrote all interactions to file: {}".format(fname))
