#!/usr/bin/env python
import argparse
import gzip
import diachrscripts_classes as dclass
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import binom



parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--fpkm-file', help='Path to a \'*_genes.fpkm_tracking\' file produced with cuffdiff.')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--use-linear-regression', help='', default='false', choices=['true','false'])

args = parser.parse_args()
out_prefix = args.out_prefix
fpkm_tracking_file = args.fpkm_file
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file
use_linear_regression = args.use_linear_regression

def categorizeDigestMaxApproach(digest):
    """
    Returns the highest expression level category for given digest.

    :param chr_name:
    :param d_sta:
    :param d_end:
    :return: Expression level category of the digest. Either 0, 1 or None.
    """
    digest_expression_categories = []

    # iterate digest and collect expression categories of TSS on digest
    for i in range(digest.get_start(), digest.get_end()):
        key = digest.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_category(key) != -1: # '-1' means no annotated TSS at this position
            digest_expression_categories.append(ref_gene_tss_map.get_coord_category(key))

    if 1 in digest_expression_categories:
        return 1 # one active TSS is enough to make the coordinate active
    elif 0 in digest_expression_categories:
        return 0
    else:
        return None # the digest has annotated TSS but FPKM is missing

def categorizeDigestPair(interaction):
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
    digest_category_1 = categorizeDigestMaxApproach(interaction.get_first_digest())
    digest_category_2 = categorizeDigestMaxApproach(interaction.get_second_digest())
    return str(digest_category_1) + "/" + str(digest_category_2)

def getInteractionTypeBinom(simple, twisted):
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

### Start Execution
###################

ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

ref_gene_tss_map.parse_cuffdiff_genes_fpkm_tracking_file(fpkm_tracking_file) # add FPKM values for genes
ref_gene_tss_map.set_expression_categories(1)

# iterate over interaction file and determine counts of pair categories
strand_simple = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_twisted = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_undirected = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_undefined = dclass.PairKeyDict(['0', '1', '-1', 'None'])

n_interaction = 0
print "[INFO] Determining pair category for each interaction..."
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    line = fp.readline()

    while line:

        # parse line representing one interaction
        values = line.split("\t")

        digest_1 = dclass.Digest(values[0], int(values[1]), int(values[2]))
        if values[3] == 'A':
            digest_1.set_active()
        digest_2 = dclass.Digest(values[4], int(values[5]), int(values[6]))
        if values[7] == 'A':
            digest_2.set_active()

        if ':' in values[8]: # regular Diachromatic file
            values2 = values[8].split(":")
            n_simple = int(values2[0])
            n_twisted = int(values2[1])
            i_type = "TBD" # will be determined using binomial P-value
        else: # from LRT script
            n_simple = values[8]
            n_twisted = values[9]
            i_type = values[13].rstrip()
        interaction = dclass.Interaction(digest_1, digest_2, n_simple, n_twisted)
        interaction.set_interaction_type(i_type)

        # restrict analysis to cis interactions
        if not(interaction.is_cis()):
            line = fp.readline()
            continue

        # restrict analysis to long range interactions
        if interaction.get_digest_distance() < 10000:
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if interaction.get_status_pair_flag() != "AA":
            line = fp.readline()
            continue

        # assign expression level category to digest using max approach
        pair_key = categorizeDigestPair(interaction)

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "S":
            strand_simple.pair_dict[pair_key] = strand_simple.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "T":
            strand_twisted.pair_dict[pair_key] = strand_twisted.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "U":
            strand_undirected.pair_dict[pair_key] = strand_undirected.pair_dict[pair_key] + 1
        elif interaction.get_interaction_type() == "NA":
            strand_undefined.pair_dict[pair_key] = strand_undefined.pair_dict[pair_key] + 1
        else:
            line = fp.readline()
            print interaction.get_interaction_type()
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T' or 'U' but was " + interaction.get_interaction_type() + "."
                            )

        n_interaction += 1
        if n_interaction%1000 == 0:
            print "\t[INFO]", n_interaction, "interactions processed..."

        line = fp.readline()

    print("...done.")

fp.close()

# print results to screen
print out_prefix
print "PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tUNDEFINED" # absolute numbers
for i in strand_simple.pair_dict.keys():
    print i + "\t" + str(strand_simple.pair_dict[i]) + "\t" + str(strand_twisted.pair_dict[i]) + "\t" + str(strand_simple.pair_dict[i] + strand_twisted.pair_dict[i]) + "\t" + str(strand_undirected.pair_dict[i]) + "\t" + str(strand_undefined.pair_dict[i])

print "PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tUNDEFINED" # relative frequencies within simple, twisted and undirected
for i in strand_simple.pair_dict.keys():
    print i + "\t" + str(float(1.0 * strand_simple.pair_dict[i] / sum(strand_simple.pair_dict.values()))) + "\t" + str(float(1.0 * strand_twisted.pair_dict[i] / sum(strand_twisted.pair_dict.values()))) + "\t" \
          + str(float((1.0 * strand_simple.pair_dict[i] + strand_twisted.pair_dict[i]) / (sum(strand_simple.pair_dict.values()) + sum(strand_twisted.pair_dict.values())))) + "\t" \
          + str(float(1.0 * strand_undirected.pair_dict[i] / sum(strand_undirected.pair_dict.values()))) + "\t" \
          + str(float(1.0 * strand_undefined.pair_dict[i] / sum(strand_undefined.pair_dict.values())))


### Graveyard
#############

# def categorizeDigestPairLinearRegression(interaction, gene_id_to_fpkm_hash):
#
#     chr_name_1 = interaction.get_first_digets().get_chromosome()
#     d_sta_1 = interaction.get_first_digets().get_start()
#     d_end_1 = interaction.get_first_digets().get_end()
#
#     chr_name_2 = interaction.get_second_digets().get_chromosome()
#     d_sta_2 = interaction.get_second_digets().get_start()
#     d_end_2 = interaction.get_second_digets().get_end()
#
#     # get FPKM values for both digest
#     FPKM_D1 = []
#     FPKM_D2 = []
#     for i in range(d_sta_1, d_end_1):
#         key = chr_name_1 + ":" + str(i)
#         if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in gene_id_to_fpkm_hash:
#             FPKM_D1.append(float(gene_id_to_fpkm_hash[tss_pos_to_gene_id[key]]))
#
#     for i in range(d_sta_2, d_end_2):
#         key = chr_name_2 + ":" + str(i)
#         if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in gene_id_to_fpkm_hash:
#             FPKM_D2.append(float(gene_id_to_fpkm_hash[tss_pos_to_gene_id[key]]))
#     xx = []
#     yy = []
#     for i in FPKM_D1:
#         for j in FPKM_D2:
#             xx = xx[:len(FPKM_D2)-1]
#             yy = yy[:len(FPKM_D2)-1]
#             xx.append(i)
#             yy.append(j)
#             x = np.array(xx).reshape((-1, 1))
#             y = np.array(yy)
#             model_1 = LinearRegression().fit(x, y)
#             r_sq = model_1.score(x, y)
#             if r_sq < 1 and 0.2 < r_sq:
#                 print "======"
#                 print FPKM_D1
#                 print FPKM_D2
#                 print len(FPKM_D1)
#                 print len(FPKM_D2)
#                 print len(xx)
#                 print len(xx)
#                 print "x:", x
#                 print "y:", y
#                 print('Coefficient of determination D1:', r_sq)
#                 xx *= 0
#                 yy *= 0
#
#     # # perform LR to determine category of D2
#     # x = np.array(FPKM_D2).reshape((-1,1))
#     # y = np.array(FPKM_D1)
#     # model_2 = LinearRegression().fit(x, y)
#     # r_sq = model_2.score(x, y)
#     #print('Coefficient of determination D2:', r_sq)
#
#     return "Foo"


# def categorizeDigest(chr_name, d_sta, d_end):
#     """
#     Input parameters are coordinates of a digest. This function checks for each digest position whether there is
#     a TSS with associated expression level category (either 0 or 1 for active/inactive or 0 to 4).
#     If all associated expression level categories are the same, the function returns the corresponding category.
#     Otherwise, if there are discordant categories, the this function returns a 'd'.
#     Finally, if there is no TSS on the digest at all, the function returns a '-1'.
#     This may happen due to slight differences of the annotation used for intercation calling, RNA-seq analysis
#     and refGene file. Ideally, the same annotation is used for all analyses.
#
#     :param chr_name: Chromosome name, e.g. chr1
#     :param d_sta: First position of a digest.
#     :param d_end: Last position of a digest.
#     :return: Concordant expression level category or 'd' if discordant categories were found or '-1' of no TSS was found.
#     """
#     current_expression_category = -1
#     for i in range(d_sta, d_end):
#         key = chr_name + ":" + str(i)
#         if key in tss_pos_to_gene_id and tss_pos_to_gene_id[key] in expression_categories:
#             if current_expression_category == -1: # first TSS on digest
#                 current_expression_category = expression_categories[tss_pos_to_gene_id[key]]
#                 continue
#             if current_expression_category != expression_categories[tss_pos_to_gene_id[key]]:
#                 return "d"
#
#     return current_expression_category