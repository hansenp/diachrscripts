#!/usr/bin/env python

from collections import defaultdict
import argparse
import gzip
import requests
import os
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
import wget

parser = argparse.ArgumentParser(description='Calculate mean distance to nearest CTCF site for sets of digests')
parser.add_argument('-i', help='input digest file')

ENSEMBL_REG_BUILD_URL = 'ftp://ftp.ensembl.org/pub/release-99/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz'

args = parser.parse_args()
infile = args.i

print("infile: %s" % infile)


class CTCF:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __overlaps(self, other_start, other_end):
        """
        Determine whether the "other" interval overlaps
        This function should only be used by 'distance'
        When we get here, we know the chromosome matches
        :param start: start position on the chromosome
        :param end: end position on the chromosome
        :return: true if other overlaps with this CTCF
        """
        if other_start > self.end:
            return False
        if other_end < self.start:
            return False
        # if we get here, we now that other_start is located upstream of self.end
        # and other_end is located downstream of self.start
        # therefore, there is some form of overlap
        return True

    def distance(self, other_chrom, other_start, other_end):
        if self.chrom != other_chrom:
            raise ValueError("Chromosomes do not match for distance calculation. Should never happen, Bug!")
        if self.__overlaps(other_start, other_end):
            return 0
        if not isinstance(other_start, int):
            print("OTHER START NOT AN INT:", other_start)
        if not isinstance(self.end, int):
            print("self.end NOT AN INT:", self.end)
        if other_start > self.end:
            return other_start - self.end
        if not isinstance(other_end, int):
            print("OTHER END NOT AN INT:", other_start)
        if not isinstance(self.start, int):
            print("self.start NOT AN INT:", self.end)
        if other_end < self.start:
            return self.start - other_end
        # we should never get here
        else:
            raise ValueError("Bad values for distance calculation")


class Digest:
    def __init__(self, chromstring, type):
        """
        Parse chr8:38379664-38385061 to 8 (no 'chr') and the two positions
        :param chromstring: e.g., chr8:38379664-38385061
        :param type:
        """
        if chromstring.startswith("chr"):
            chromstring = chromstring[3:]
        A = chromstring.split(":")
        if len(A) != 2:
            raise TypeError("Bad chrstr (A)", chromstring)
        self.chrom = A[0]
        B = A[1].split("-")
        if len(B) != 2:
            raise TypeError("Bad chrstr (B)", chromstring)
        self.start = int(B[0])
        self.end = int(B[1])
        self.type = type


def extract_ctcf_from_ensembl(url):
    local_name = 'homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz'
    local_path = os.path.join('tmp_data', local_name)
    ctcf_binding_sites = []
    if os.path.exists(local_path):
        print("regulatory build file previously downloaded!")
    else:
        print("Downloading regulatory build file")
        wget.download(url, local_path)
    # When we get here, the ensembl file is available
    with gzip.open(local_path, 'rt') as g:
        for line in g:
            # print(line)
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                print("[ERROR] Bad line", line)
                continue
            chrom = fields[0]
            categ = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            if "CTCF" in categ:
                ctcf = CTCF(chrom, start, end)
                ctcf_binding_sites.append(ctcf)
    return ctcf_binding_sites


def process_digests():
    digest_list = []
    n_chromosome_Y = 0
    with gzip.open(infile, 'rt') as g:
        for line in g:
            # print(line)
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                raise ValueError("Bad digest line", line)
            loc = fields[0]
            type = fields[2]
            chroms = loc.split(';')
            if loc.startswith('chrY'):
                # Note that the regulatory build does not have chromosome Y
                n_chromosome_Y += 1
                continue
            if len(chroms) != 2:
                raise ValueError("Bad chroms line=", line)
            digest_a = Digest(chroms[0], type)
            digest_list.append(digest_a)
            digest_b = Digest(chroms[1], type)
            digest_list.append(digest_b)
    print(
        "Skipped %d chromosome Y digest pairs (Ensembl regulatory build does not include chromosome Y)" % n_chromosome_Y)
    return digest_list


ctcf_binding_sites = extract_ctcf_from_ensembl(ENSEMBL_REG_BUILD_URL)
ctcf_dict = defaultdict(list)
for ctcf in ctcf_binding_sites:
    ctcf_dict[ctcf.chrom].append(ctcf)
print("Extracted %d CTCF binding sites" % len(ctcf_binding_sites))
for k, v in ctcf_dict.items():
    print("Chrom {}: {} CTCF sites".format(k, len(v)))
digest_list = process_digests()
print("Extracted %d digests" % len(digest_list))
digest_dict = defaultdict(list)
for d in digest_list:
    digest_dict[d.type].append(d)
for k, v in digest_dict.items():
    print("Category {}: {} digests".format(k, len(v)))


# When we get here, all of the data are ingested and we want to calculate the distance. We are doing this
# in a naive way, but it should work

def calculate_distances(ctcf_d, digest_list):
    """
    calculate the distance of each digest to the nearest CTCF
    :param ctcf_d: complete dictionary of CTCF sites according to chromosome
    :param digest_list: one list of digests, e.g. DIII
    :return: list of distances
    """
    min_dist = []
    n = 0
    BIG_NUMBER = 1000000000
    LIMIT = 10000  # for testing
    for d in digest_list:
        chr = d.chrom
        start = d.start
        end = d.end
        ctcf_list = ctcf_d[chr]
        if ctcf_list is None:
            raise ValueError("Could not get CTCF list for chrom {}".format(chr))
        distance = BIG_NUMBER
        for ctcf in ctcf_list:
            try:
                dist = ctcf.distance(chr, start, end)
            except Error:
                continue
            if dist < distance:
                distance = dist
        if distance == BIG_NUMBER:
            print("Dist is BUG")
            print("Digest=chr{}:{}-{}".format(d.chrom, d.start, d.end))
        min_dist.append(distance)
        n += 1
        if n > LIMIT:
           break
        if n % 1000 == 0:
            print("Processed %d digests" % n)
    return np.array(min_dist)


for k, v in digest_dict.items():
    print("Category {}: {} digests".format(k, len(v)))
    min_dist = calculate_distances(ctcf_dict, v)
    print("\t{} +/- {}".format(np.mean(min_dist), np.std(min_dist)))
    N = len(min_dist)
    oneslist = [1 for i in min_dist if i > 0]
    n_overlap = sum(oneslist)
    print("\t{}/{} ({:.2f}%) digests had overlapping CTCF sites".format(n_overlap, N, 100.0*n_overlap/N))
