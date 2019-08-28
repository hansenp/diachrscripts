import argparse
import gzip
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file


### Define auxiliary functions
##############################

def get_digest_strand_tag(digest, ref_gene_tss_map):
    strand = -1
    for i in range(digest.get_start(), digest.get_end()): # iterate digest from left to right
        key = digest.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_strand(key) == -1: # no TSS at this position
            continue
        elif ref_gene_tss_map.get_coord_strand(key) == 'd': # TSS for '-' and '+' at this same position (rare)
            return 'd'
        else:
            if strand == -1: # first TSS on digest
                strand = ref_gene_tss_map.get_coord_strand(key) # must be '-' or '+'
            elif strand != ref_gene_tss_map.get_coord_strand(key):
                return 'd'
    return strand

def get_interaction_strand_pair_key(interaction, ref_gene_tss_map):
    d1_strand = get_digest_strand_tag(interaction.get_first_digest(), ref_gene_tss_map)
    d2_strand = get_digest_strand_tag(interaction.get_second_digest(), ref_gene_tss_map)
    return str(d1_strand) + "/" + str(d2_strand)


### Start execution
###################

# prepare variables and data structures
ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

strand_simple = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_twisted = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_undirected = dclass.PairKeyDict(['+', '-', 'd', '-1'])
strand_indefinable = dclass.PairKeyDict(['+', '-', 'd', '-1'])

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

        if n_interaction_total%10000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(10000)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if interaction.get_digest_status_pair_flag() != "AA":
            n_non_promoter_promoter_interaction +=1
            line = fp.readline()
            continue

        # assign expression level category to digest using max approach
        pair_key = get_interaction_strand_pair_key(interaction, ref_gene_tss_map)

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
            strand_indefinable.pair_dict[pair_key] = strand_indefinable.pair_dict[pair_key] + 1 # less than 5 read pairs
            n_indefinable_interaction += 1
        else:
            line = fp.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

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
print( out_prefix)
print("PAIR\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tINDEFINABLE\tSIMPLE\tTWISTED\tSIMPLE+SIMPLE\tUNDIRECTED\tINDEFINABLE")
for key in strand_simple.pair_dict.keys():

    # absolute frequencies
    print("\'" + key + "\'" + "\t" + str(strand_simple.pair_dict[key]) + "\t" + str(strand_twisted.pair_dict[key]) + "\t" + str(strand_simple.pair_dict[key] + strand_twisted.pair_dict[key]) + "\t" + str(strand_undirected.pair_dict[key]) + "\t" + str(strand_indefinable.pair_dict[key]))

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
