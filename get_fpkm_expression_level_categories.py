#!/usr/bin/env python
import argparse
import gzip
import sys
import numpy as np
from scipy.stats import binom

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.')
parser.add_argument('--fpkm-file', help='Path to a \'*_genes.fpkm_tracking\' file produced with cuffdiff.')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be the same that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--categorization-model', help='Choose \'two\' for inactive/active categorization and \'five\' for categorization with five expression levels.')

args = parser.parse_args()
out_prefix = args.out_prefix
fpkm_tracking_file = args.fpkm_file
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file
categorization_model = args.categorization_model

def init_pair_hashs(mode):
    """
    :param mode: Either 'two' for simple active/inactive categorization or 'five' for categorization with five expression levels.
    :return: Nothing. Only the global hashes for digest category pair counting are initialized.
    """
    if mode == "two":
        for i in ['d','0','1','-1']:
            for j in ['d','0','1','-1']:
                key = i + "/" + j
                PAIR_hash_simple[key] = 0
                PAIR_hash_twisted[key] = 0
                PAIR_hash_undirected[key] = 0
    if mode == "five":
        for i in ['d','0','1','2','3','4','-1']:
            for j in ['d','0','1','2','3','4','-1']:
                key = i + "/" + j
                PAIR_hash_simple[key] = 0
                PAIR_hash_twisted[key] = 0
                PAIR_hash_undirected[key] = 0

def categorizeDigest(chr_name, d_sta, d_end):
    """
    Input parameters are coordinates of a digest. This function checks for each digest position whether there is
    a TSS with associated expression level category (either 0 or 1 for active/inactive or 0 to 4).
    If all associated expression level categories are the same, the function returns the corresponding category.
    Otherwise, if there are discordant categories, the this function returns a 'd'.
    Finally, if there is no TSS on the digest at all, the function returns a '-1'.
    This may happen due to slight differences of the annotation used for intercation calling, RNA-seq analysis
    and refGene file. Ideally, the same annotation is used for all analyses.

    :param chr_name: Chromosome name, e.g. chr1
    :param d_sta: First position of a digest.
    :param d_end: Last position of a digest.
    :return: Concordant expression level category or 'd' if discordant categories were found or '-1' of no TSS was found.
    """
    current_expression_category = -1
    for i in range(d_sta, d_end):
        key = chr_name + ":" + str(i)
        if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in expression_categories:
            if current_expression_category == -1: # first TSS on digest
                current_expression_category = expression_categories[tss_pos_to_gene_id[key]]
                continue
            if current_expression_category != expression_categories[tss_pos_to_gene_id[key]]:
                return "d"

    return current_expression_category

def categorizeDigestPair(chr_name_1, d_sta_1, d_end_1, chr_name_2, d_sta_2, d_end_2):
    """
    This function uses the function 'categorizeDigest' in order to determine the expression level categories of two
    interacting digests and assembles the corresponding pair key, e.g. 1/1 if both digests have category 1.

    :param chr_name_1: Chromosome name of first digest.
    :param d_sta_1: First position of a first digest.
    :param d_end_1: Last position of first digest.
    :param chr_name_2: Chromosome name of second digest.
    :param d_sta_2: First position of a second digest.
    :param d_end_2: Last position of second digest.
    :return:
    """
    digest_category_1 = categorizeDigest(chr_name_1, d_sta_1, d_end_1)
    digest_category_2 = categorizeDigest(chr_name_2, d_sta_2, d_end_2)
    return str(digest_category_1) + "/" + str(digest_category_2)

def getInteractionType(simple,twisted):
    """
    This function determines whether a given interaction is significantly directed, i.e. simple or twisted,
    and returns a 'S' for simple and 'T' for twisted. Otherwise, the function returns an 'U' for undirected
    or a 'NA', if the sum of twisted and simple is smaller than five because interactions with less read pairs
    cannot be significant using the binomial distribution.

    :param simple: Number of simple read pairs.
    :param twisted: Number of Twisted read pairs.
    :return: 'S' (simple),'T' (twisted), 'U' (undirected), 'NA' (not available)
    """
    if simple + twisted < 5:
        return "NA"

    if simple < twisted:
        p_value = 1 - binom.cdf(twisted - 1, simple + twisted, 0.5)
        if p_value <= 0.05:
            return "T"
        else:
            return "U"
    else:
        p_value = 1 - binom.cdf(simple - 1, simple + twisted, 0.5)
        if p_value <= 0.05:
            return "S"
        else:
            return "U"

def parse_cuffdiff_genes_fpkm_tracking_file(fpkm_tracking_file):
    """
    This function parses a *_genes.fpkm_tracking file produced by 'cuffdiff' and determines the number of genes with
    zero FPKM as well as the upper limits of the first three quartiles used for categorization into five expression
    levels (or inctive/active categories).

    :param fpkm_tracking_file: TSV file with FPKM values for genes
    :return: Number of genes with zero FPKM and upper limits of the first three quartiles.
    """
    n_zero_fpkm = 0

    # read non-zero FPKM vaules to array
    fpkm_array = []
    with open(fpkm_tracking_file) as fp:
        line = fp.readline()

        while line:
            values = line.split("\t")

            if values[0] == "tracking_id": # skip first line
                line = fp.readline()
                continue
            else:
                fpkm = float(values[9])
                if 0 < fpkm:
                    fpkm_array.append(fpkm)
                else:
                    n_zero_fpkm += 1

            line = fp.readline()
    fp.close()

    # determine quartiles
    upper_first_q = np.quantile(fpkm_array, .25)
    upper_second_q = np.quantile(fpkm_array, .50)
    upper_third_q = np.quantile(fpkm_array, .75)

    return n_zero_fpkm, upper_first_q, upper_second_q, upper_third_q

def get_expression_level_category_hash(fpkm_tracking_file, model, upper_first_q, upper_second_q, upper_third_q):
    """
    This function assigns each gene an expression level category based on FPKM values. There are two models for
    categorization: For the model 'two', genes are either inactive (0) or active (1). For the model 'five',
    genes with zero FPKM belong to category 0, and the remaining genes to category 1 to 4 depending on FPKM quartiles.

    :param fpkm_tracking_file: fpkm_tracking_file: TSV file with FPKM values for genes
    :param upper_first_q: Upper limit of first quartile.
    :param upper_second_q: Upper limit of second quartile.
    :param upper_third_q: Upper limit of third quartile.
    :return: Hash with gene IDs as keys and category numbers as values.
    """
    expression_categories = {}

    with open(fpkm_tracking_file) as fp:
        line = fp.readline()

        while line:
            values = line.split("\t")

            if values[0] == "tracking_id":
                line = fp.readline()
                continue

            gene_id = values[3]
            fpkm = float(values[9])

            if model == "two": # only inactive/active
                if fpkm == 0:
                    expression_categories[gene_id] = 0
                elif fpkm < upper_first_q:
                    expression_categories[gene_id] = 0
                elif fpkm < upper_second_q:
                    expression_categories[gene_id] = 1
                elif fpkm < upper_third_q:
                    expression_categories[gene_id] = 1
                else:
                    expression_categories[gene_id] = 1
            elif model == "five":
                if fpkm == 0:
                    expression_categories[gene_id] = 0
                elif fpkm < upper_first_q:
                    expression_categories[gene_id] = 1
                elif fpkm < upper_second_q:
                    expression_categories[gene_id] = 2
                elif fpkm < upper_third_q:
                    expression_categories[gene_id] = 3
                else:
                    expression_categories[gene_id] = 4
            else:
                raise Exception('[FATAL] Invalid categorization model. Should be either \'two\' or \'five\' but was {}.', model)

            line = fp.readline()

    fp.close()
    return expression_categories

def parse_refGene_file(ref_gene_file):

    tss_pos_to_gene_id = {}

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

            line = fp.readline()

    fp.close()

    return tss_pos_to_gene_id

### Start Execution
###################

# determine genes with zero FPKM and quartiles for categorization into expression levels
n_zero_fpkm, upper_first_q, upper_second_q, upper_third_q = parse_cuffdiff_genes_fpkm_tracking_file(fpkm_tracking_file)

print "[INFO] There were", n_zero_fpkm, "genes with zero FPKM."
print "[INFO] Upper FPKM limit of the first quartile:", upper_first_q
print "[INFO] Upper FPKM limit of the second quartile:", upper_second_q
print "[INFO] Upper FPKM limit of the third quartile:", upper_third_q

# init hash map, key=<gene_id, e.g. NR_157147>, value=<category, i.e. 0 and 1 or 0 to 4>
expression_categories = get_expression_level_category_hash(fpkm_tracking_file, "two", upper_first_q, upper_second_q, upper_third_q)

# init hash map, key=<coordinates, e.g. chr1:12345>, value=<gene_id, e.g. NR_157147>
tss_pos_to_gene_id = parse_refGene_file(ref_gene_file)

# iterate over interaction file and determine counts of pair categories
PAIR_hash_simple = {} # Keys: digest pair categories; Values: Corresponding counts
PAIR_hash_twisted = {}
PAIR_hash_undirected = {}
init_pair_hashs("two")

with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
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

        if len(values2) == 2: # regular Diachromatic file
            cnt_simple = int(values2[0])
            cnt_twisted = int(values2[1])
            interaction_type = getInteractionType(cnt_simple, cnt_twisted)
        else: # file from LRT script
            interaction_type = values[13]

        pair_key = categorizeDigestPair(chr_name_1,d_sta_1,d_end_1,chr_name_2, d_sta_2, d_end_2)

        if interaction_type == "NA":
            line = fp.readline()
            continue
        elif interaction_type == "S\n" or interaction_type == "S":
            PAIR_hash_simple[pair_key] = PAIR_hash_simple[pair_key] + 1
        elif interaction_type == "T\n" or interaction_type == "T":
            PAIR_hash_twisted[pair_key] = PAIR_hash_twisted[pair_key] + 1
        elif interaction_type == "U\n" or interaction_type == "U":
            PAIR_hash_undirected[pair_key] = PAIR_hash_undirected[pair_key] + 1
        else:
            line = fp.readline()
            print interaction_type
            exit("Interaction type is neither NA, S, T or U. This should never happen.")

        line = fp.readline()

fp.close()

print out_prefix
# print results to screen
print "SIMPLE\tTWISTED\tUNDIRECTED" # absolute numbers
for i in PAIR_hash_simple:
    print i + "\t" + str(PAIR_hash_simple[i]) + "\t" + str(PAIR_hash_twisted[i]) + "\t" + str(PAIR_hash_undirected[i])

print "SIMPLE\tTWISTED\tUNDIRECTED" # relative frequencies within simple, twisted and undirected
for i in PAIR_hash_simple:
    print i + "\t" + str(float(1.0*PAIR_hash_simple[i]/sum(PAIR_hash_simple.values()))) + "\t" + str(float(1.0*PAIR_hash_twisted[i]/sum(PAIR_hash_twisted.values()))) + "\t" + str(float(1.0*PAIR_hash_undirected[i]/sum(PAIR_hash_undirected.values())))
