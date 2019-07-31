#!/usr/bin/env python
import sys
import numpy as np
from scipy.stats import binom

name = sys.argv[1] # genes.fpkm_tracking file of cuffdiff
fpkm_tracking_file = sys.argv[2] # genes.fpkm_tracking file of cuffdiff
ref_gene_file = sys.argv[3] # UCSC refGene file (must be the same that was used to create the digest map for Diachromatic)
diachromatic_interaction_file = sys.argv[4] # UCSC refGene file (must be the same that was used to create the digest map for Diachromatic)

PAIR_hash_simple = {}
PAIR_hash_twisted = {}
PAIR_hash_undirected = {}

def getDigestCategory(chr_name, d_sta, d_end):
    for i in range(d_sta, d_end):
        key=chr_name + ":" + str(i)
        if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in expression_categories:
            print tss_pos_to_gene_id[key], "\t", expression_categories[tss_pos_to_gene_id[key]], "\t", tss_pos_to_strand[key]

def getDigestCategoryStrandPairs(chr_name_1, d_sta_1, d_end_1, chr_name_2, d_sta_2, d_end_2, direction):
    D1_Array = []
    for i in range(d_sta_1, d_end_1):
        key=chr_name_1 + ":" + str(i)
        if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in expression_categories:
            D1_Array.append(str(expression_categories[tss_pos_to_gene_id[key]]) + tss_pos_to_strand[key]) # collect TSS and expression for first digest
    for i in range(d_sta_2, d_end_2):
        key = chr_name_2 + ":" + str(i)
        if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in expression_categories:
            for j in D1_Array:
                pair_key = j + "/" + str(expression_categories[tss_pos_to_gene_id[key]]) + tss_pos_to_strand[key]
                if direction == "simple":
                    if pair_key in PAIR_hash_simple:
                        PAIR_hash_simple[pair_key] = PAIR_hash_simple[pair_key] + 1
                    else:
                        PAIR_hash_simple[pair_key] = 1
                elif direction == "twisted":
                    if pair_key in PAIR_hash_twisted:
                        PAIR_hash_twisted[pair_key] = PAIR_hash_twisted[pair_key] + 1
                    else:
                        PAIR_hash_twisted[pair_key] = 1
                else:
                    if pair_key in PAIR_hash_undirected:
                        PAIR_hash_undirected[pair_key] = PAIR_hash_undirected[pair_key] + 1
                    else:
                        PAIR_hash_undirected[pair_key] = 1

# read non-zero FPKM vaules to array
fpkm_array = []
with open(fpkm_tracking_file) as fp:
    line = fp.readline()

    while line:
        values = line.split("\t")

        if values[0] == "tracking_id":
            line = fp.readline()
            continue
        else:
            fpkm = float(values[9])
            if 0 < fpkm:
                fpkm_array.append(fpkm)

        line = fp.readline()

fp.close()

# determine and report quartiles
print "Quartiles for", name, ":", np.quantile(fpkm_array, (.25, .5,.75))
upper_first_q = np.quantile(fpkm_array,.25)
upper_second_q = np.quantile(fpkm_array,.50)
upper_third_q = np.quantile(fpkm_array,.75)

# init hash map, key=<gene_id, e.g. NR_157147>, value=<category, i.e. 0 to 4>
expression_categories = {}

# reread genes.fpkm_tracking file and fill hash hap
with open(fpkm_tracking_file) as fp:
    line = fp.readline()

    while line:
        values = line.split("\t")

        if values[0] == "tracking_id":
            line = fp.readline()
            continue

        gene_id = values[3]
        fpkm = float(values[9])
        if fpkm == 0:
            expression_categories[gene_id] = 0
        elif fpkm < upper_first_q:
            expression_categories[gene_id] = 0
        elif fpkm < upper_second_q:
            expression_categories[gene_id] = 0
        elif fpkm < upper_third_q:
            expression_categories[gene_id] = 1
        else:
            expression_categories[gene_id] = 1

        line = fp.readline()

fp.close()

# init hash map, key=<coordinates, e.g. chr1:12345>, value=<gene_id, e.g. NR_157147>
tss_pos_to_gene_id = {}
tss_pos_to_strand = {}

# read refGene.txt.gz file to hash map
with open(ref_gene_file) as fp:
    line = fp.readline()

    while line:
        values = line.split("\t")

        gene_id = values[1]
        chr_name = values[2]
        tss_strand = values[3]
        if tss_strand == "+":
            tss_pos = values[4]
        else:
            tss_pos = values[5]

        tss_pos_to_gene_id[chr_name + ":" + tss_pos] = gene_id
        tss_pos_to_strand[chr_name + ":" + tss_pos] = tss_strand

        line = fp.readline()

fp.close()

# iterate over interaction file
with open(diachromatic_interaction_file) as fp:
    line = fp.readline()

    while line:
        values = line.split("\t")

        chr_name_1 = values[0]
        d_sta_1 = int(values[1])
        d_end_1 = int(values[2])
        d_state_1 = values[3]

        chr_name_2 = values[4]

        d_sta_2 = int(values[5])
        d_end_2 = int(values[6])

        d_state_2 = values[7]

        # restrict analysis to cis interactions (why?)
        if chr_name_1 != chr_name_2:
            line = fp.readline()
            continue

        # restrict analysis to long range interactions
        if d_sta_2 - d_end_1 < 10000:
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if d_state_1 != "A" or d_state_2 != "A":
            line = fp.readline()
            continue

        # calculate binomial P-value
        values2 = values[8].split(":")
        cnt_simple = int(values2[0])
        cnt_twisted = int(values2[1])
        cnt_total = cnt_simple + cnt_twisted

        if cnt_total < 5:
            line = fp.readline()
            continue

        if cnt_simple < cnt_twisted:
            p_value = 1 - binom.cdf(cnt_twisted-1, cnt_total, 0.5)
        else:
            p_value = 1 - binom.cdf(cnt_simple-1, cnt_total, 0.5)

        if p_value <= 0.05:  # significantly directed interactions


            if cnt_twisted < cnt_simple: # simple
                getDigestCategoryStrandPairs(chr_name_1,d_sta_1,d_end_1,chr_name_2, d_sta_2, d_end_2, "simple")
            else:  # twisted
                getDigestCategoryStrandPairs(chr_name_1, d_sta_1, d_end_1, chr_name_2, d_sta_2, d_end_2, "twisted")
        else:  # undirected (not significantly directed)
            getDigestCategoryStrandPairs(chr_name_1, d_sta_1, d_end_1, chr_name_2, d_sta_2, d_end_2, "undirected")

        line = fp.readline()

fp.close()

print "SIMPLE"
for i in PAIR_hash_simple:
    print i + "\t" + str(PAIR_hash_simple[i])

print "TWISTED"
for i in PAIR_hash_twisted:
    print i + "\t" + str(PAIR_hash_twisted[i])

print "UNDIRECTED"
for i in PAIR_hash_undirected:
    print i + "\t" + str(PAIR_hash_undirected[i])
