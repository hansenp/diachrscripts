#!/usr/bin/env python

"""
This script takes a Diachromatic interaction file and a 'refGene.txt.gz' file and creates for each interaction a line in
enhanced interaction format with the binomial P-value for the counts of simple and twisted read pairs, interaction
distances and, if the associated digests contain TSS, TSS coordinates and gene symbols as well as strand information.

This is one line in enhanced interaction format as an example:

chr2:112534779-112543248;chr2:112577153-112587985	33905	NA	LOC105373562,POLR1B;CHCHD5	221:130	AA	14.20	d/+	chr2:112541915:+,chr2:112542212:-,chr2:112542036:+;chr2:112584437:+,chr2:112584609:+,chr2:112584854:+

Each individual field of an enhanced interaction line is explained below.

------------------------------------------------------------------------------------------------------------------------

FIELD 1: 'chr2:112534779-112543248;chr2:112577153-112587985'

Coordinates of the interacting digests separated
separated by a semicolon. The first comes before the second digest in sequential order.

------------------------------------------------------------------------------------------------------------------------

FIELD 2: '33905'

Distance between interacting digests measured as the end position of the first an start position of
the second digest.

------------------------------------------------------------------------------------------------------------------------

FIELD 3: 'NA'

Wildcard for interaction category tag.

------------------------------------------------------------------------------------------------------------------------

FIELD 4: 'LOC105373562,POLR1B;CHCHD5'

Two comma separated lists of gene symbols separated by a semicolon. The symbols
before and after the semicolon correspond to the TSS on the first and second digest.

------------------------------------------------------------------------------------------------------------------------

FIELD 5: '221:130'

Number of simple and twisted read pair counts separated by a colon.

------------------------------------------------------------------------------------------------------------------------

FIELD 6: 'AA'

Enrichment pair tag. Indicates which digest of the interaction was selected for target enrichment,
whereby 'A' means 'active' (enriched) and 'I' 'inactive' (not enriched). There are four possible tags:

   1. 'AA' - Both digests selected
   2. 'AI' - First digests selected
   3. 'IA' - Second digests selected
   4. 'II' - No digests selected

------------------------------------------------------------------------------------------------------------------------

FIELD 7: '14.20'

Negative logarithm of the binomial P-value for orientation of interactions.

------------------------------------------------------------------------------------------------------------------------

FIELD 8: 'd/+'

Strand pair tag. There are three strand tags for digests:

   1. '-' - Digest contains one or more TSS on the minus strand only
   2. '+' - Digest contains one or more TSS on the plus strand only
   3. 'd' - Digest contains TSS on the minus and plus strand (discordant)

At the level of interactions, this results in nine possible strand pair tags.

------------------------------------------------------------------------------------------------------------------------

FIELD 9: 'chr2:112541915:+,chr2:112542212:-,chr2:112542036:+;chr2:112584437:+,chr2:112584609:+,chr2:112584854:+'

Two comma separated lists of TSS coordinates separated by a semicolon.

------------------------------------------------------------------------------------------------------------------------
"""


import argparse
import gzip
import math
import diachrscripts_toolkit as dclass


### Parse command line
######################

parser = argparse.ArgumentParser(description='Create enhanced interaction file from Diachromatic interaction and refGene.txt.gz file.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--ref-gene-file', help='UCSC refGene file (must be gzipped and the same version that was used to create the digest map for Diachromatic).')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
args = parser.parse_args()

out_prefix = args.out_prefix
ref_gene_file = args.ref_gene_file
diachromatic_interaction_file = args.interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + diachromatic_interaction_file)
print("\t[INFO] refGene file: " + ref_gene_file)


### Define auxiliary functions
##############################

def get_gene_symbols_of_interacting_digests(interaction, ref_gene_tss_map):
    """
    This function takes an interaction object and a map containing information about genes, i.e. gene symbols as well as
    coordinates and strands of TSS. For each digest of the interaction, the associated annotation data is retrieved from
    the map and combined into arrays.

    :param interaction: Interaction object of the class 'Interaction' defined in 'diachrscripts_toolkit.py'
    :param ref_gene_tss_map: Object of the class 'TSSCoordinateMap' defined in 'diachrscripts_toolkit.py'

    :return: gene_symbols_d1: Array containing all symbols of genes on the first digest
             gene_symbols_d2: Array containing all symbols of genes on the second digest
             gene_strands_d1: Array containing all strand tags of TSS ('-', '+' or 'd') on the first digest
             gene_strands_d2: Array containing all strand tags of TSS ('-', '+' or 'd') on the second digest
             tss_list_d1: Array containing all coordinates of TSS ('<chr>:<pos>:<strand>') on the first digest
             tss_list_d2: Array containing all coordinates of TSS ('<chr>:<pos>:<strand>') on the second digest
    """

    gene_symbols_d1 = []
    gene_strands_d1 = []
    tss_list_d1 = []
    d1 = interaction.get_first_digest()
    for i in range(d1.get_start(), d1.get_end() + 1):  # iterate digest from left to right
        key = d1.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_symbols(key) == -1:  # no TSS at this position
            continue
        else:
            gene_symbols_d1.append(ref_gene_tss_map.get_coord_symbols(key))
            gene_strands_d1.append(list(ref_gene_tss_map.get_coord_strand(key)))
            tss_list_d1.append(key + ":" + ref_gene_tss_map.get_coord_strand(key))

    gene_symbols_d2 = []
    gene_strands_d2 = []
    tss_list_d2 = []
    d2 = interaction.get_second_digest()
    for i in range(d2.get_start(), d2.get_end() + 1):  # iterate digest from left to right
        key = d2.get_chromosome() + ":" + str(i)
        if ref_gene_tss_map.get_coord_symbols(key) == -1:  # no TSS at this position
            continue
        else:
            gene_symbols_d2.append(ref_gene_tss_map.get_coord_symbols(key))
            gene_strands_d2.append(list(ref_gene_tss_map.get_coord_strand(key)))
            tss_list_d2.append(key + ":" + ref_gene_tss_map.get_coord_strand(key))

    return gene_symbols_d1, gene_symbols_d2, gene_strands_d1, gene_strands_d2, tss_list_d1, tss_list_d2


### Prepare variables, data structures and streams for output files
###################################################################

file_name_original = out_prefix + "_original_interactions_with_genesymbols.tsv"
f_output_original = open(file_name_original, 'wt')

ref_gene_tss_map = dclass.TSSCoordinateMap(ref_gene_file, "refGene") # parse refGene file with TSS
ref_gene_tss_map.analyze_coordinates_and_print_report() # collect counts and print report


### Start execution
###################

print("[INFO] Adding P-value, gene symbols and TSS coordinates to each interaction in " + diachromatic_interaction_file + " ...")

print("\t[INFO] Iterating Diachromatic interaction file ...")

with gzip.open(diachromatic_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Create interaction object from Diachromatic interaction line
        interaction = dclass.Interaction(line)

        # Calculate logarithm of binomial P-value
        interaction.set_logsf_binomial_p_value()
        logsf_binomial_p_value = "{:.2f}".format(-interaction.get_binomial_p_value())

        # Fetch gene symbols and TSS coordinates from refGene.txt.gz
        d1_symbols, d2_symbols, d1_strands, d2_strands, tss_list_d1, tss_list_d2 = get_gene_symbols_of_interacting_digests(interaction, ref_gene_tss_map)

        # Combine gene symbols into comma separated lists separated by a semicolon
        symbols_d12 = ",".join(set(sum(d1_symbols, [] ))) + ";" + ",".join(set(sum(d2_symbols, [] )))

        # Combine strands of TSS into comma separated lists separated by a semicolon
        strands_d12 = ",".join(set(sum(d1_strands, []))) + ";" + ",".join(set(sum(d2_strands, [])))

        # Combine TSS coordinates into comma separated lists separated by a semicolon
        tss_d12 = ",".join(set(tss_list_d1)) + ";" + ",".join(set(tss_list_d2))

        # Determine strand category tag ('-', '+' or 'd') of the first digest
        sd1='-1'
        if ['d'] in d1_strands:
            sd1='d'
        elif ['+'] in d1_strands and ['-'] in d1_strands:
            sd1 = 'd'
        elif ['+'] in d1_strands:
            sd1 = '+'
        elif ['-'] in d1_strands:
            sd1 = '-'

        # Determine strand category tag ('-', '+' or 'd') of the second digest
        sd2='-1'
        if ['d'] in d2_strands:
            sd2='d'
        elif ['+'] in d2_strands and ['-'] in d2_strands:
            sd2 = 'd'
        elif ['+'] in d2_strands:
            sd2 = '+'
        elif ['-'] in d2_strands:
            sd2 = '-'

        # Combine strand tags into strand pair tag
        strand_pair_tag = sd1 + "/" + sd2

        # Combine simple and twisted read pairs into one colon separated string
        simple_twisted_counts = str(interaction.n_simple) + ":" + str(interaction.n_twisted)

        # Get distance between interacting digests
        d_dist = str(interaction.get_digest_distance())

        # Set wildcard for interaction category
        itype = "NA"

        # Write interaction to enhanced interaction file
        f_output_original.write(
            interaction.get_coord_string() + "\t" +
            d_dist  + "\t" +
            itype + "\t" +
            symbols_d12 + "\t" +
            simple_twisted_counts + "\t" +
            interaction.get_digest_pair_flag_original_order() + "\t" +
            logsf_binomial_p_value + "\t" +
            strand_pair_tag + "\t" +
            tss_d12
            + "\n")

        line = fp.readline()

print("\t[INFO] ... done.")

fp.close()
f_output_original.close()

print("[INFO] Number of processed interactions: " + str(n_interaction_total))
