#!/usr/bin/env python
import argparse
import gzip
import diachrscripts_toolkit as dclass


### Parse command line
######################

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


### Define auxiliary functions
##############################

def get_digest_expression_tag(digest, ref_gene_tss_map):
    """
    This function traverses the input digest from left to right and, for each position, queries a help structure
    that allows assign the digest a expression tag.

    :param digest: A digest that is part of an interaction.
    :param ref_gene_tss_map: A help structure that combines information about TSS and expression.
    :return: Expression level tag of the digest. Either 0 (inactive), 1 (active) or None (?).
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

def get_interaction_expression_tag_pair_key(interaction, ref_gene_tss_map):
    digest_category_1 = get_digest_expression_tag(interaction.get_first_digest(), ref_gene_tss_map)
    digest_category_2 = get_digest_expression_tag(interaction.get_second_digest(), ref_gene_tss_map)
    return str(digest_category_1) + "/" + str(digest_category_2)


### Start execution
###################

ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

ref_gene_tss_map.parse_cuffdiff_genes_fpkm_tracking_file(fpkm_tracking_file) # add FPKM values for genes
ref_gene_tss_map.set_expression_categories(2) # set expression categories depending on FPKM quartiles

# iterate over interaction file and determine counts of pair categories
strand_simple = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_twisted = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_undirected = dclass.PairKeyDict(['0', '1', '-1', 'None'])
strand_indefinable = dclass.PairKeyDict(['0', '1', '-1', 'None'])

n_interaction_total = 0
n_trans_short_range_interaction = 0
n_non_promoter_promoter_interaction = 0
n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_indefinable_interaction = 0

# iterate interactions
print("[INFO] Determining pair category for each interaction in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:

    line = fp.readline()

    while line:

        if n_interaction_total%1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)
        interaction.set_interaction_type("TBD", 0.003)

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(10000)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if interaction.get_digest_status_pair_flag() != "AA":
            n_non_promoter_promoter_interaction += 1
            line = fp.readline()
            continue

        # assign expression level category to digest using max approach
        pair_key = get_interaction_expression_tag_pair_key(interaction, ref_gene_tss_map)

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "S":
            strand_simple.pair_dict[pair_key] = strand_simple.pair_dict[pair_key] + 1
            n_simple_interaction += 1
        elif interaction.get_interaction_type() == "T":
            strand_twisted.pair_dict[pair_key] = strand_twisted.pair_dict[pair_key] + 1
            n_twisted_interaction += 1
        elif interaction.get_interaction_type() == "U":
            strand_undirected.pair_dict[pair_key] = strand_undirected.pair_dict[pair_key] + 1
            n_undirected_interaction += 1
        elif interaction.get_interaction_type() == "NA":
            strand_indefinable.pair_dict[pair_key] = strand_indefinable.pair_dict[pair_key] + 1
            n_indefinable_interaction += 1
        else:
            line = fp.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T' or 'U' but was " + interaction.get_interaction_type() + ".")

        line = fp.readline()

    print("... done.")

fp.close()


### Print results to screen
###########################

# sum up simple and twisted counts
simple_twisted_dict = {}
for key in strand_simple.pair_dict.keys():
    simple_twisted_dict[key] = strand_simple.pair_dict[key] + strand_twisted.pair_dict[key]

# print absolute frequencies, calculate and print relative frequencies
print(out_prefix)
print("PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tINDEFINABLE\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tINDEFINABLE")
for key in strand_simple.pair_dict.keys():

    # absolute frequencies
    print("\'" + key + "\'" + "\t" + str(strand_simple.pair_dict[key]) + "\t" + str(strand_twisted.pair_dict[key]) + "\t" + str(strand_simple.pair_dict[key] + strand_twisted.pair_dict[key]) + "\t" + str(strand_undirected.pair_dict[key]) + "\t" + str(strand_indefinable.pair_dict[key])),

    # relative frequencies
    fraction_simple = dclass.get_fraction(strand_simple.pair_dict[key], strand_simple.pair_dict.values())
    fraction_twisted = dclass.get_fraction(strand_twisted.pair_dict[key], strand_twisted.pair_dict.values())
    fraction_directed = dclass.get_fraction(simple_twisted_dict[key], simple_twisted_dict.values())
    fraction_undirected = dclass.get_fraction(strand_undirected.pair_dict[key], strand_undirected.pair_dict.values())
    fraction_indefinable = dclass.get_fraction(strand_indefinable.pair_dict[key], strand_indefinable.pair_dict.values())
    print("\t" + fraction_simple + "\t" + fraction_twisted + "\t" + str(fraction_directed) + "\t" + str(fraction_undirected) + "\t" + str(fraction_indefinable))

print("Total number of interactions: " + str(n_interaction_total))
print("Number of trans and short range interactions: " + str(n_trans_short_range_interaction) + " (discarded)")
print("Number of non promoter-promoter interactions: " + str(n_non_promoter_promoter_interaction) + " (discarded)")
print("Number of directed simple interactions: " + str(n_simple_interaction))
print("Number of directed twisted interactions: " + str(n_twisted_interaction))
print("Number of undirected interactions: " + str(n_undirected_interaction))
print("Number of indefinable interactions: " + str(n_indefinable_interaction))


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