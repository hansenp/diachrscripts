import argparse
import gzip
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--p-value-cutoff', help='P-value cutoff used for categorization of interactions.')

args = parser.parse_args()
out_prefix = args.out_prefix
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file
p_value_cutoff = float(args.p_value_cutoff)

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + diachromatic_interaction_file)
print("\t[INFO] refGene file: " + ref_gene_file)
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))


### Define auxiliary functions
##############################

def get_gene_symbols_of_interacting_digests(interaction, ref_gene_tss_map):

    gene_symbols_d1 = []
    d1 = interaction.get_first_digest()
    for i in range(d1.get_start(), d1.get_end()):  # iterate digest from left to right
        key = d1.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_symbols(key) == -1:  # no TSS at this position
            continue
        else:
            gene_symbols_d1.append(ref_gene_tss_map.get_coord_symbols(key))

    gene_symbols_d2 = []
    d2 = interaction.get_second_digest()
    for i in range(d2.get_start(), d2.get_end()):  # iterate digest from left to right
        key = d2.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_symbols(key) == -1:  # no TSS at this position
            continue
        else:
            gene_symbols_d2.append(ref_gene_tss_map.get_coord_symbols(key))

    return gene_symbols_d1, gene_symbols_d2


### Start execution
###################

file_name_original = out_prefix + "_original_interactions_with_genesymbols.tsv"
f_output_original = open(file_name_original, 'wt')

# prepare variables and data structures
ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report

n_interaction_total = 0
n_trans_short_range_interaction = 0
n_non_promoter_promoter_interaction = 0
n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_indefinable_interaction = 0

# Given the P-value cutoff, find the smallest n that yields a significant P-value
n_indefinable_cutoff = dclass.find_indefinable_n(p_value_cutoff)

# Determine distribution of n for directed interactions
min_digest_dist = 10000
status_pair_flag = "AA"
n_dict = dclass.get_n_dict(diachromatic_interaction_file, status_pair_flag, min_digest_dist, p_value_cutoff)

# iterate interactions
print("[INFO] Determining pair category for each interaction in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, 'rt') as fp:

    line = fp.readline()

    while line:

        if n_interaction_total%10000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(min_digest_dist)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        # restrict analysis to interactions between targeted promoters
        if interaction.get_digest_status_pair_flag() != status_pair_flag:
            n_non_promoter_promoter_interaction +=1
            line = fp.readline()
            continue

        # set the type of interaction based on P-value ('S', 'T', 'U', 'NA')
        if interaction.get_interaction_type() == "TBD":
            interaction.set_interaction_type("TBD", p_value_cutoff, n_indefinable_cutoff)

        # assign expression level category to digest using max approach
        d1_symbols, d2_symbols = get_gene_symbols_of_interacting_digests(interaction, ref_gene_tss_map)

        symbols_d12 = ",".join(set(sum(d1_symbols, [] ))) + ";" + ",".join(set(sum(d2_symbols, [] )))

        simple_twisted_counts = str(interaction.n_simple) + ":" + str(interaction.n_twisted)

        n_total = interaction.n_simple + interaction.n_twisted

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "NA":
            n_indefinable_interaction += 1
            itype = "NA"
        elif interaction.get_interaction_type() == "U":
            n_undirected_interaction += 1
            if n_total in n_dict and 0 < n_dict[n_total]:
                itype = "UR"
                n_dict[n_total] = n_dict[n_total] - 1
            else:
                itype = "U"
        elif interaction.get_interaction_type() == "S":
            n_simple_interaction += 1
            itype = "S"
        elif interaction.get_interaction_type() == "T":
            n_twisted_interaction += 1
            itype = "T"
        else:
            line = fp.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        line = line.rstrip()

        f_output_original.write(interaction.get_coord_string() + "\t" + str(interaction.get_digest_distance())  + "\t" + itype + "\t" + symbols_d12 + "\t" + simple_twisted_counts + "\n")

        line = fp.readline()

    print("... done.")

fp.close()
f_output_original.close()