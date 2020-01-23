"""
This script performs a number of comparisons between directed and undirected interactions and associated
digests and promoters.

It builds upon the output of the script 'get_gene_symbols_of_interactions_script.py' that contains predefined
sets of directed and undirected reference interactions (column 3).

The analysis is restricted to the following interaction categories:

   1. Directed interactions between promoters: Simple (S) or twisted (T) and in addition enrichment status pair tag AA
   2. Undirected reference interactions between promoters (URAA)

Directed and undirected reference interactions were selected in way so as they are comparable with regard to:

   1. Signal strength measured as the number of read pairs
   2. Enrichment status of interacting digests (AA)
   3. Distance between interacting digests



# Overview
##########

This script performs performs comparisons of interactions, digests and associated promoters between directed (D) and
undirected reference (UR) interactions. For this purpose, comparative sets for directed and undirected interactions need to be
defined at the level of digests and interactions. The comparative sets can be defined in two different ways:

1. Compare digests that occur exclusively either for directed or undirected interactions, i.e. remove all ambiguous digests that
are involved in directed and undirected interactions at the same time. We refer to the remaining digest as
'exclusive digests' (ED). At the level of interactions, keep only interactions with two exclusive digests. We refer to
these interaction as 'exclusive interactions' (EI). This yields four comparative sets:

   i. Exclusive directed digests EDD
   ii. Exclusive directed digests EUD
   iii. Exclusive directed interactions EDI
   iv. Exclusive directed interactions EUI

2. Remove ambiguous digests for undirected interactions only and leave digests from directed interactions unfiltered.
This yields four comparative sets:

   i. Directed digests DD
   ii. Exclusive directed digests EUD
   iii. Directed interactions DI
   iv. Exclusive directed interactions EUI

Use the option --comparative-sets with 'DWOU-VS-UWOD'  for the first and 'D-VS-UWOD' for the second definition of
comparative sets.



# Input file
############

The input file was prepared using the python script 'get_gene_symbols.py' and contains
detailed information about one interaction in each line. This is one line as an example:

chr2:112534779-112543248;chr2:112577153-112587985	33905	S	LOC105373562,POLR1B;CHCHD5	221:130	AA	14.20	d/+	chr2:112541915:+,chr2:112542212:-,chr2:112542036:+;chr2:112584437:+,chr2:112584609:+,chr2:112584854:+

Field 1: 'chr2:112534779-112543248;chr2:112577153-112587985' - Coordinates of the interacting digests separated
separated by a semicolon. The first comes before the second digest in sequential order.

Field 2: '33905' - Distance between interactiong digests measured as the end position of the first an start position
of the second digest.

Field 3: 'S' - Tag for interaction category. 'S' means directed simple interactions. Beyond that, there are
'T' for directed twisted, 'U' undirected, 'NA' for indefinable and the undirected reference interactions
'URII', 'URAI' and 'URAA' for interacions with no, one or two enriched digests. The categorization is done by
the script 'get_gene_symbols.py' based on a P-value threshold defined using the empirical FDR procedure. For this
script, we use typically 'S', 'T' and 'URAA' because we are dealing with promoter-promoter interactions.

Field 4: 'LOC105373562,POLR1B;CHCHD5' - Two comma separated lists of gene symbols separated by a semicolon. The symbols
before and after the semicolon correspond to the TSS on the first and second digest of the interaction.

Field 5: '221:130' - Number of simple and twisted read pair counts separated by a colon.

Field 6: 'AA' - Enrichment pair tag. Indicates which digest of the interaction was selected for target enrichment,
whereby 'A' means 'active' (enriched) and 'I' 'inactive' (not enriched). There are four possible tags:
'AA' (both digests selected), 'AI' (first digests selected), 'IA' (second digests selected), 'II' (no digests selected).

Field 7: '14.20' - Negative decadic logarithm of the binomial P-value for orientation of interactions.

Field 8: 'd/+' - Strand pair tag. There are three tags for the individual strands: '-', '+' and 'd' separted by a
forward slash. The first tag corresponds to the first and the second tag to the second digest of the interactions.
A digest is assigned a '-' tag, if it has one or more TSS on the reverse strand and no TSS on the forward strand.
The same applies to digest that are assigned a '+'. Digests wit TSS on the forward and reverse strand are assigned an
'd' (stands for discordant).

Field 9: 'chr2:112541915:+,chr2:112542212:-,chr2:112542036:+;chr2:112584437:+,chr2:112584609:+,chr2:112584854:+' - Two
comma separated lists of TSS (chromosome:coordinate:strand) separated by a semicolon.

# Workflow
##########

Essentially, the script follows a two pass approach:

   1. Pass through all interactions in order to identify unique exclusive digests that do not interact with digests
   from the other group, i.e. remove digests from directed interactions, if they additionally interact with a digest
   from an undirected reference interaction, and remove digests from undirected reference interactions, if they
   additionally interact with a digest from a directed interaction.

   2. Pass through all interactions a second time in order to identify TSS that are on unique exclusive digests.


# Calculated statistics
#######################

In the course of that, a variety of digest and promoter features is determined that are listed below:

1. Connectivity factor: One digest can be involved in more than one interaction. During the first pass we use a set to
store the coordinates of digests from directed interactions and undirected reference interactions. In order to
calculate the connectivity factor for directed and undirected references interactions, we divide the cardinalities of the
corresponding sets by twice the number of corresponding interactions and subtract this number from 1.

2. Mean number of interaction partners of unique exclusive digests: These numbers are determined during the second
pass.

3. Base frequencies:

4. Digest sizes:

5. Interaction sizes:

6. Strand pair tag frequencies:


# Output files
##############

This script generates the following files:

# Central file with most of the calculated statistics in one tab separated line

<OUT_PREFIX>_interaction_and_digest_statistics.tab

# Base and strand pair tag frequencies

<OUT_PREFIX>_base_frequencies.tab
<OUT_PREFIX>_strand_pair_tag_frequencies.tab

# Gene symbols associated with digests for GO analysis

<OUT_PREFIX>_directed_digest_symbols.tab
<OUT_PREFIX>_undirected_digest_symbols.tab


# Exhaustive collections for interaction partners per digests and digest and interaction sizes

<OUT_PREFIX>_directed_digests_interaction_partner_per_digest_distribution.tab
<OUT_PREFIX>_undirected_digests_interaction_partner_per_digest_distribution.tab

<OUT_PREFIX>_directed_digests_size_distribution.tab
<OUT_PREFIX>_undirected_digests_size_distribution.tab

<OUT_PREFIX>_directed_interaction_size_distribution.tab
<OUT_PREFIX>_undirected_interaction_size_distribution.tab


# Generated BED and FASTA files for calculation of sequence statistics and downstream analysis:

<OUT_PREFIX>_directed_digests.bed                                                    # Digest regions (keep for others)
<OUT_PREFIX>_directed_digests.fasta                                                  # Digest sequences for statistics (REMOVED after calculation)
<OUT_PREFIX>_directed_digests_tss_redundant.bed                                      # Redundant promoter regions (keep for others)
<OUT_PREFIX>_directed_digests_tss_redundant_merged.bed                               # Merged promoter regions (keep for others)
<OUT_PREFIX>_directed_digests_tss_redundant_merged.fasta                             # Merged promoter sequences for statistics (REMOVED after calculation)
<OUT_PREFIX>_directed_digests_tss_redundant_merged_center_trimmed.bed                # Merged and center-trimmed promoter regions (keep for others)
<OUT_PREFIX>_directed_digests_tss_redundant_merged_center_trimmed_masked.fasta       # Merged and center-trimmed promoter repeat masked sequences for motif analysis with DREME
<OUT_PREFIX>_directed_digests_tss_redundant_minus.fasta                              # Redundant promoter regions on '-' strand for motif analysis with CentriMo
<OUT_PREFIX>_directed_digests_tss_redundant_plus.fasta                               # Redundant promoter regions on '+' strand for motif analysis with CentriMo

<OUT_PREFIX>_undirected_digests.bed                                                  # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests.fasta                                                # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant.bed                                    # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_merged.bed                             # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_merged.fasta                           # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_merged_center_trimmed.bed              # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_merged_center_trimmed_masked.fasta     # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_minus.fasta                            # see comment for directed counterpart
<OUT_PREFIX>_undirected_digests_tss_redundant_plus.fasta                             # see comment for directed counterpart
"""

import argparse
import gzip
import os
import numpy


### Parse command line
######################

parser = argparse.ArgumentParser(description='Extract regions for motif analysis from file with interactions and gene symbols.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-gs-file', help='Interaction file with gene symbols created with \'get_gene_symbols_of_interactions_script.py\'.')
parser.add_argument('--chrom-info-file', help='File with chromosome sizes <CHROMOSOME NAME>\\t<CHROMOSOME SIZE>.')
parser.add_argument('--up-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--down-dist', help='Number of bases upstream of TSS.', default=1000)
parser.add_argument('--allowed-enrichment-pair-tags', help='Set of allowed enrichment pair tags for digests (\'AA\', \'AI\',\'IA\',\'II\', ...) separated by \';\'.', default="AA")
parser.add_argument('--allowed-strand-pair-tags', help='Set of allowed strand pair tags for digests (\'-/-\', \'+/+\', \'+/-\',\'-/+\',\'-/d\', ...) separated by \';\'.', default="-/-;+/+")
parser.add_argument('--allowed-interaction-categories-directed', help='Set categories of directed interactions (\'S\' and \'T\') separated by \';\'.', default="S;T")
parser.add_argument('--allowed-interaction-categories-undirected', help='Set categories of directed interactions (\'NA\', \'U\', \'URAA\', \'URAI\' and \'URII\') separated by \';\'.', default="URAA")
parser.add_argument('--comparative-sets', help='Sets of digests/interactions that will be compared. Two options directed vs. undirected without directed (D-VS-UWOD) or directed without undirected vs. undirected without directed (DWOU-VS-UWOD).', default="D-VS-UWOD")
parser.add_argument('--bedtools-path', help='Path to bedtools executable.', default="bedtools")
parser.add_argument('--genome-fasta-path', help='Path to genome FASTA file.', default="/Users/hansep/data/hg38/hg38.fa")

args = parser.parse_args()
out_prefix = args.out_prefix
interaction_gs_file = args.interaction_gs_file # central file that contains a lot of information about interactions
chrom_info_file = args.chrom_info_file # needed to avoid promoter coordinates outside chromosome regions
up_dist = int(args.up_dist) # used to determine start coordinates of promoter regions
down_dist = int(args.down_dist) # used to determine end coordinates of promoter regions
allowed_enrichment_pair_tags = args.allowed_enrichment_pair_tags # only interactions with these enrichment pair tags will be used
allowed_strand_pair_tags = args.allowed_strand_pair_tags # only interactions with these strand pair tags will be used
allowed_interaction_categories_directed = args.allowed_interaction_categories_directed # only directed interactions with these categories will be used
allowed_interaction_categories_undirected = args.allowed_interaction_categories_undirected # only undirected with these categories will be used
comparative_sets = args.comparative_sets
bedtools_path = args.bedtools_path
genome_fasta_path = args.genome_fasta_path

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] --interaction-gs-file: " + interaction_gs_file)
print("\t[INFO] --chrom-info-file: " + chrom_info_file)
print("\t[INFO] --up-dist: " + str(up_dist))
print("\t[INFO] --down-dist: " + str(down_dist))
print("\t[INFO] --enrichment-strand-pair-tags: " + allowed_enrichment_pair_tags)
print("\t[INFO] --allowed-strand-pair-tags: " + allowed_strand_pair_tags)
print("\t[INFO] --allowed-interaction-categories-directed: " + allowed_interaction_categories_directed)
print("\t[INFO] --allowed-interaction-categories-undirected: " + allowed_interaction_categories_undirected)
print("\t[INFO] --comparative-sets: " + comparative_sets)
print("\t[INFO] --bedtools-path: " + bedtools_path)
print("\t[INFO] --genome-fasta-path: " + genome_fasta_path)

# convert allowed pair tag strings into lists
allowed_strand_pair_tags = allowed_strand_pair_tags.split(";")
allowed_enrichment_pair_tags = str(allowed_enrichment_pair_tags).split(";")
allowed_interaction_categories_directed = str(allowed_interaction_categories_directed).split(";")
allowed_interaction_categories_undirected = str(allowed_interaction_categories_undirected).split(";")


### Define auxiliary functions
##############################

def convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name):
    """
    This function uses the external tool 'bedtools' in order to extract sequences whose coordinates are provided in
    BED formatted file.

    :param bedtools_path: Path to external command line tool 'bedtools'
    :param genome_fasta_path: Path to genome FASTA file containing one sequence for eaach chromosome
    :param bed_file_name: BED file with coordinates of the sequences that are going to be extracted
    :return: Filename of the created FASTA file
    """
    bed_file_name_base = str(bed_file_name).split(".b")[0]
    sys_cmd = bedtools_path + ' getfasta -name -fi ' + genome_fasta_path + ' -bed  ' + bed_file_name + ' > ' + bed_file_name_base + '.fasta'
    os.system(sys_cmd)
    return bed_file_name_base + '.fasta'


def merge_bed(bedtools_path, bed_file_name):
    """
    This function uses the external tool 'bedtools' in order to merge the regions of a BED file.

    :param bedtools_path: Path to external command line tool 'bedtools'
    :param bed_file_name: BED file to be merged
    :return: Filename of merged BED file
    """
    bed_file_name_base = str(bed_file_name).split(".b")[0]
    sys_cmd = bedtools_path + ' sort -i ' + bed_file_name + ' > ' + bed_file_name_base + '_sorted.bed'
    os.system(sys_cmd)
    sys_cmd = bedtools_path + ' merge -i ' + bed_file_name_base + '_sorted.bed' + ' > ' + bed_file_name_base + '_merged.bed'
    os.system(sys_cmd)
    sys_cmd = 'rm ' + bed_file_name_base + '_sorted.bed'
    os.system(sys_cmd)
    return bed_file_name_base + '_merged.bed'


def center_trim_bed(bed_file_name, up_size, down_size):
    """
    Takes a BED file and trims regions that are longer than 'up_size + down_size' coordinates of the trimmed regions
    are 'region_center - up_size' and 'region_center - down_size'.
    If a region is shorter that 'up_size + down_size', a warning will be thrown.

    :param bed_file_name: BED file that is going to be trimmed
    :param up_size: Distance from the center in upstream direction
    :param down_size: Distance from the center in downstream direction
    :return:
        Total number of regions in input BED file,
        number of regions that are longer than 'up_size + down_size',
        size of largest region in input BED file and
        name of the created BED file with trimmed regions
    """
    bed_file_name_base = str(bed_file_name).split(".b")[0]
    bed_file_name_center_trimmed = bed_file_name_base + '_center_trimmed.bed'
    bed_file_stream_center_trimmed = open(bed_file_name_center_trimmed, 'wt')

    target_region_size = up_size + down_size
    total_region_num = 0                       # total number of regions in input BED file
    trimmed_region_num = 0                     # number of trimmed regions that are longer than 'up_size + down_size'
    region_max_size = 0                        # size of largest region in input BED file
    total_region_size = 0                      # total size of regions in input BED file
    total_trimmed_region_size = 0              # total size of trimmed region

    with open(bed_file_name, 'rt') as fp:
        line = fp.readline()
        while line:
            total_region_num += 1

            chr = line.split("\t")[0]
            sta = int(line.split("\t")[1])
            end = int(line.split("\t")[2])
            size = end - sta

            total_region_size = total_region_size + size

            if region_max_size < size:
                region_max_size = size

            if target_region_size < size: # trim
                trimmed_region_num += 1
                center = sta + int(size / 2)
                sta = center - up_size
                end = center + down_size
                total_trimmed_region_size = total_trimmed_region_size + (size - target_region_size)

            bed_file_stream_center_trimmed.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(total_region_num) + "\n")

            if size < target_region_size:
                print("[WARNING] Region is smaller than \'up_dist + down_dist\'!")

            line = fp.readline()

    bed_file_stream_center_trimmed.close()

    avg_region_size = int(total_region_size / total_region_num)                                 # Average size of input regions
    trimmed_region_num_frac = "{:.2f}".format(100*trimmed_region_num / total_region_num)                # fraction of input regions that were trimmed
    total_trimmed_region_size_frac = "{:.2f}".format(100*total_trimmed_region_size/total_region_size)   # overall fraction of bases that were trimmed (removed)

    return target_region_size, total_region_num, trimmed_region_num_frac, avg_region_size, region_max_size, total_trimmed_region_size_frac, bed_file_name_center_trimmed


def mask_repeats(fasta_file_name):
    """
    This function takes a FASTA file and replaces all occurrences of small letters (a,c,g and t) with N.

    :param fasta_file_name: FASTA formatted file
    :return: Filename of the created repeat masked FASTA file with suffix '_masked.fasta'
    """
    fasta_file_name_base = str(fasta_file_name).split(".f")[0]
    sys_cmd = 'awk \'{if($1 !~ /^>/){gsub(/a|c|g|t/,"N",$1)};print}\' ' + fasta_file_name + ' > ' + fasta_file_name_base + '_masked.fasta'
    os.system(sys_cmd)
    return fasta_file_name_base + '_masked.fasta'


def get_base_frequencies(fasta_file_name):
    """
    This function takes a FASTA file and counts the overall number of occurrences of a, c, g, t, A, C, G, T and N.
    It calculates the repeat content as measured by the fraction of small letters (soft masked repeats).
    Furthermore, it calculates the GC content within repeat regions, outside repeat regions and the overall GC content.

    :param fasta_file_name: FASTA file with soft masked repeats
    :return: header_line: header line for a table that contains all results including absolute counts
             value_line: values that correspond to the header line
             repeat_content: Overall proportion of small letters (N regions included)
             gc_content_repeat: Proportion of g and c within repeat regions
             gc_content_non_repeat: Proportion of G and C outside repeat regions
             gc_content_total: Overall proportion of g, c, G and C (N regions excluded)
    """

    # Iterate FASTA file and count letter occurrences
    a_abs = 0
    c_abs = 0
    g_abs = 0
    t_abs = 0
    A_abs = 0
    C_abs = 0
    G_abs = 0
    T_abs = 0
    N_abs = 0
    with open(fasta_file_name, 'rt') as fp:
        line = fp.readline()
        while line:
            if line.count('>') == 0:
                line = line.rstrip()
                a_abs = a_abs + line.count('a')
                c_abs = c_abs + line.count('c')
                g_abs = g_abs + line.count('g')
                t_abs = t_abs + line.count('t')
                A_abs = A_abs + line.count('A')
                C_abs = C_abs + line.count('C')
                G_abs = G_abs + line.count('G')
                T_abs = T_abs + line.count('T')
                N_abs = N_abs + line.count('N')
            line = fp.readline()
    fp.close()

    # Calculate repeat content (fraction of lower case)
    acgt_abs = a_abs + c_abs + g_abs + t_abs
    acgtACGTN_abs = a_abs + c_abs + g_abs + t_abs + A_abs + C_abs + G_abs + T_abs + N_abs
    repeat_content = "{:.2f}".format(100*acgt_abs/acgtACGTN_abs)

    # Determine GC content within repeat regions
    cg_abs = c_abs + g_abs
    if 0 < acgt_abs:
        gc_content_repeat = "{:.2f}".format(100 * cg_abs / acgt_abs)
    else:
        gc_content_repeat = "{:.2f}".format(0.0)

    # Determine GC content for non repeat regions
    CG_abs = C_abs + G_abs
    ACGT_abs = A_abs + C_abs + G_abs + T_abs
    gc_content_non_repeat = "{:.2f}".format(100 * CG_abs / ACGT_abs)

    # Determine overall GC content
    cgCG_abs = cg_abs + CG_abs
    acgtACGT_abs = a_abs + c_abs + g_abs + t_abs + A_abs + C_abs + G_abs + T_abs
    gc_content_total = "{:.2f}".format(100 * cgCG_abs / acgtACGT_abs)

    # Prepare header line
    header_line = "fasta_file_name" + "\t" + \
                  "a_occ_abs" + "\t" + \
                  "c_occ_abs" + "\t" + \
                  "g_occ_abs" + "\t" + \
                  "t_occ_abs" + "\t" + \
                  "A_occ_abs" + "\t" + \
                  "C_occ_abs" + "\t" + \
                  "G_occ_abs" + "\t" + \
                  "T_occ_abs" + "\t" + \
                  "N_occ_abs" + "\t" + \
                  "repeat_content" + "\t" + \
                  "gc_content_repeat" + "\t" + \
                  "gc_content_non_repeat" + "\t" + \
                  "gc_content_total"

    # Prepare corresponding line with values
    value_line = fasta_file_name + "\t" + \
                 str(a_abs) + "\t" + \
                 str(c_abs) + "\t" + \
                 str(g_abs) + "\t" + \
                 str(t_abs) + "\t" + \
                 str(A_abs) + "\t" + \
                 str(C_abs) + "\t" + \
                 str(G_abs) + "\t" + \
                 str(T_abs) + "\t" + \
                 str(N_abs) + "\t" + \
                 str(repeat_content) + "\t" + \
                 str(gc_content_repeat) + "\t" + \
                 str(gc_content_non_repeat) + "\t" + \
                 str(gc_content_total)

    return header_line, value_line, repeat_content, gc_content_repeat, gc_content_non_repeat, gc_content_total


def determine_digest_sizes_and_write_to_file(tab_file_name_digests_size_distribution, digest_set):
    """
    This function takes a file name and a set of digest coordinates and writes all digest sizes to a file.
    Furthermore it calculates and returns quantiles (Q1, Q2 and Q3) and mean digest sizes.

    :param tab_file_name_digests_size_distribution: File name for a text file that will contain one digest size per line
    :param digest_set: Set of digest coordinates given as tab separated string
    :return: q1: 0.25 quantile of size distribution
             q2: 0.50 quantile of size distribution (median)
             q3: 0.75 quantile of size distribution
             mean: Average digest size
    """

    # Open stream for output file
    tab_stream_name_digests_size_distribution = open(tab_file_name_digests_size_distribution, 'wt')

    # Iterate digest set, add sizes to array and write sizes to file
    digest_size_array = []
    for d_coord in digest_set:
        sta = d_coord.split("\t")[1]
        end = d_coord.split("\t")[2]
        size = int(end) - int(sta)
        digest_size_array.append(size)
        tab_stream_name_digests_size_distribution.write(str(size) + "\n")

    # Close stream for output file
    tab_stream_name_digests_size_distribution.close()

    # Calculate quantiles and mean
    q1 = int(numpy.quantile(digest_size_array, 0.25))
    q2 = int(numpy.quantile(digest_size_array, 0.5))
    q3 = int(numpy.quantile(digest_size_array, 0.75))
    mean = int(numpy.mean(digest_size_array))

    return q1, q2, q3, mean


def write_number_of_interaction_partners_per_digest(tab_file_name_digests_interaction_partner_per_digest_distribution, interaction_partners_per_digest_dict):
    """
    This function takes a dictionary that contains the number of interaction partners of digests and prints one number
    per line to a file.

    :param tab_file_name_digests_interaction_partner_per_digest_distribution: Name of the file that is going to be created
    :param interaction_partners_per_digest_dict: Dictionary with digest coordinates as keys and number of interaction partners as values
    :return: mean_number_of_interaction_partners
    """

    # Open stream for output file
    tab_stream_name_digests_interaction_partner_per_digest_distribution = open(tab_file_name_digests_interaction_partner_per_digest_distribution, 'wt')

    # Create array with all numbers of interaction partners
    interaction_partner_number_array = []

    # Iterate dictionary
    for inter_partner_num in interaction_partners_per_digest_dict.values():
        tab_stream_name_digests_interaction_partner_per_digest_distribution.write(str(inter_partner_num) + "\n")
        interaction_partner_number_array.append(inter_partner_num)

    # Close stream for output file
    tab_stream_name_digests_interaction_partner_per_digest_distribution.close()

    # Calculate average number of interaction partner per digest
    mean_number_of_interaction_partners = "{:.2f}".format(numpy.mean(interaction_partner_number_array))

    return mean_number_of_interaction_partners


def get_unique_exclusive_gene_symbols_and_write_to_file(digests_symbols_set_a, digests_symbols_set_b, tab_file_name_unique_exclusive_digest_symbols, symbol_search_pattern):
    """
    This function takes two sets of gene symbols A and B and determines the set A\B, which is written to a text file
    that contains one gene symbol in each line. In addition, the number of gene symbols in A\B that contain the
    substring 'symbol_search_pattern' is determined.

    :param tab_file_name_unique_exclusive_digest_symbols: Name for the file with gene symbols that is going to be created
    :param digests_symbols_set_a: Set of gene symbols A
    :param digests_symbols_set_b: Another set of gene symbols B
    :param symbol_search_pattern: Count gene symbols containing this substring (We typically use 'ZNF')
    :return: symbols_abs: Absolute number of symbols in A\B
             symbols_search_pattern_abs: Absolute number of symbols in A\B containing the search pattern 'symbol_search_pattern'
             symbols_search_patter_rel: Relative number of symbols in A\B containing the search pattern 'symbol_search_pattern'
    """

    # Remove symbols from A that are also in B
    a_wo_b_digests_symbols_set = digests_symbols_set_a.difference(digests_symbols_set_b)

    # Open stream for output file
    tab_stream_name_unique_exclusive_digest_symbols = open(tab_file_name_unique_exclusive_digest_symbols, 'wt')

    # Iterate and count remaining symbols in A and write to file
    symbols_abs = 0
    symbols_search_pattern_abs = 0
    for symbol in a_wo_b_digests_symbols_set:
        symbols_abs += 1
        tab_stream_name_unique_exclusive_digest_symbols.write(symbol + "\n")
        if symbol_search_pattern in symbol:
            symbols_search_pattern_abs += 1

    # Close stream for output file
    tab_stream_name_unique_exclusive_digest_symbols.close()

    # Calculate proportion of gene symbols containing the search pattern
    symbols_search_patter_rel = "{:.2f}".format(100 * symbols_search_pattern_abs / symbols_abs)

    return symbols_abs, symbols_search_pattern_abs, symbols_search_patter_rel


def parse_interaction_line_with_gene_symbols(interaction_line_with_gene_symbols):
    """
    This function takes a tab separated line with coordinates and gene symbols that correspond to an interaction and
    parses it into individual fields that are relevant for the analyses performed in this script.

    :param interaction_line_with_gene_symbols: Line from a file that was created using the script 'get_gene_symbols_interactions'.
    :return: chr_a: Chromosome of the first digest of the interaction (e.g. chr1)
             sta_a: Start coordinate of the first digest of the interaction (e.g. 123456)
             end_a: End coordinate of the first digest of the interaction (e.g. 234567)
             syms_a: Comma separated list of gene symbols associated with TSS on the first digest (e.g. HIST1H1E,HIST1H2BD)
             tsss_a: Comma separated list of TSS coordinates on the first digest (e.g. chr6:26158121:+,chr6:26156331:+)
             chr_b, sta_b, end_b, syms_b, tsss_b: Same as for the first digest but for the second digest
             enrichment_pair_tag: Two letter tag indicating the enrichment status of the two digests (e.g.  AA, AI, II)
             strand_pair_tag: Two symbol tag separated by '/' indicating the strands of TSS on first and second digest (e.g. -/-, +/+, -/+, +/-, -/d, ...)
             interaction_category: Interaction category with respect to directionality (S, T, URAA, URAI, ...)
    """

    # Split line into individual fields
    field = interaction_line_with_gene_symbols.split("\t")

    # Split string for digest pair associated with the interaction
    coordinate_pair = field[0].split(";")
    coordinate_a = coordinate_pair[0]
    coordinate_b = coordinate_pair[1]

    # Extract digest coordinates
    chr_a = coordinate_a.split(":")[0]
    chr_b = coordinate_b.split(":")[0]

    sta_a = int(coordinate_a.split(":")[1].split("-")[0])
    sta_b = int(coordinate_b.split(":")[1].split("-")[0])

    end_a = int(coordinate_a.split(":")[1].split("-")[1])
    end_b = int(coordinate_b.split(":")[1].split("-")[1])

    # Split string for pair of comma separated lists of gene symbols
    symbols_pair = field[3].split(";")
    syms_a = symbols_pair[0]
    syms_b = symbols_pair[1]

    # Split string for pair of comma separated lists of TSS
    tsss_pair = field[8].split(";")
    tsss_a = tsss_pair[0]
    tsss_b = tsss_pair[1]

    # Extract letter pair tag for enrichment status of digests
    enrichment_pair_tag = field[5]

    # Extract symbol pair tag for TSS orientations on associated digests
    strand_pair_tag = field[7]

    # Extract category of interaction with respect to directionality
    interaction_category = field[2]

    return chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category


def write_interaction_sizes_to_file(interaction_sizes_array, tab_file_name_interaction_size_distribution):
    """
    This function takes an array with distances between interacting digests and prints them to file, on distance per
    line. Furthermore, it calculates and returns quantiles and mean.

    :param interaction_sizes_array: Array containing interaction sizes
    :param tab_file_name_interaction_size_distribution: Name of the file that is going to be created containing one interaction size per line
    :return: q1: 0.25 quantile of size distribution
             q2: 0.50 quantile of size distribution (median)
             q3: 0.75 quantile of size distribution
             mean: Average interaction size
    """

    # Open stream for output file
    tab_file_stream_interaction_size_distribution = open(tab_file_name_interaction_size_distribution, 'wt')

    # Iterate array with sizes and write to file
    for size in interaction_sizes_array:
        tab_file_stream_interaction_size_distribution.write(str(size) + "\n")

    # Close stream for output file
    tab_file_stream_interaction_size_distribution.close()

    # Calculate quantiles and mean
    q1_interaction_sizes = int(numpy.quantile(interaction_sizes_array, 0.25))
    q2_interaction_sizes = int(numpy.quantile(interaction_sizes_array, 0.5))
    q3_interaction_sizes = int(numpy.quantile(interaction_sizes_array, 0.75))
    mean_interaction_sizes = int(numpy.mean(interaction_sizes_array))

    return q1_interaction_sizes, q2_interaction_sizes, q3_interaction_sizes, mean_interaction_sizes


### Prepare output BED files with coordinates of digests and associated promoters
#################################################################################

# create two BED files for regions of directed and undirected interactions
bed_file_name_directed_digests = out_prefix + "_directed_digests.bed"
bed_stream_name_directed_digests = open(bed_file_name_directed_digests, 'wt')
bed_file_name_undirected_digests = out_prefix + "_undirected_digests.bed"
bed_stream_name_undirected_digests = open(bed_file_name_undirected_digests, 'wt')

# create two BED files for promoters on digests of directed iteractions, one for promoters on the plus and another one for promoters on the minus strand
bed_file_name_directed_digests_tss_minus = out_prefix + "_directed_digests_tss_redundant_minus.bed"
bed_stream_name_directed_digests_tss_minus = open(bed_file_name_directed_digests_tss_minus, 'wt')
bed_file_name_directed_digests_tss_plus = out_prefix + "_directed_digests_tss_redundant_plus.bed"
bed_stream_name_directed_digests_tss_plus = open(bed_file_name_directed_digests_tss_plus, 'wt')

# create two BED files for promoters on digests of undirected iteractions, one for promoters on the plus and another one for promoters on the minus strand
bed_file_name_undirected_digests_tss_minus = out_prefix + "_undirected_digests_tss_redundant_minus.bed"
bed_stream_name_undirected_digests_tss_minus = open(bed_file_name_undirected_digests_tss_minus, 'wt')
bed_file_name_undirected_digests_tss_plus = out_prefix + "_undirected_digests_tss_redundant_plus.bed"
bed_stream_name_undirected_digests_tss_plus = open(bed_file_name_undirected_digests_tss_plus, 'wt')

bed_file_name_directed_digests_tss = out_prefix + "_directed_digests_tss_redundant.bed"
bed_file_name_undirected_digests_tss = out_prefix + "_undirected_digests_tss_redundant.bed"


### Prepare variables and data structures
#########################################

# Read chromosome sizes from file to hash that can be used to avoid promoter coordinates outside chromosomes
chrom_sizes_dict = {}
with open(chrom_info_file, 'rt') as fp:
    line = fp.readline()
    while line:
        chr_name = line.split("\t")[0]
        chr_size = line.split("\t")[1]
        chrom_sizes_dict[chr_name] = int(chr_size)
        line = fp.readline()


# Sets that will store tab separated digest coordinates
directed_digests_set = set()
undirected_digests_set = set()

# Variables to count interactions
directed_interaction_num = 0
undirected_interaction_num = 0

# Variable to count digests that have no TSS (should never happen, if interactions are filtered for 'AA')
first_digest_without_tss_num = 0
second_digest_without_tss_num = 0

# Set that will store gene symbols
directed_digests_symbols_set = set()
undirected_digests_symbols_set = set()

# Dictionaries that will be used to count the number of interaction partners per digest
interaction_partners_per_digest_directed_dict = {}
interaction_partners_per_digest_undirected_dict = {}

# Dictionaries to count the occurences of strand pair tags ('-/-', '-/+', '+/-', ...)
strand_pair_tag_list = ["-/-", "+/+", "-/+", "+/-", "-/d", "d/-", "+/d", "d/+", "d/d"]
strand_pair_tag_directed_dict = {}
strand_pair_tag_undirected_dict = {}

### First pass: Determine unique digests that do not interact with digests from the other interaction set
#########################################################################################################

print("[INFO] First pass: Determining unique exclusive digests ...")

print("\t[INFO] Iterating interactions with gene symbols ...")

with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse line from interaction file with gene symbols
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category = \
            parse_interaction_line_with_gene_symbols(line.rstrip())

        # Skip lines/interactions with the wrong enrichment pair tag (typically we use only AA)
        if enrichment_pair_tag not in allowed_enrichment_pair_tags:
            line = fp.readline()
            continue

        # Skrip line with wrong strand symbol pair tag (typically we use no restriction, i.e. 'All')
        if strand_pair_tag not in allowed_strand_pair_tags and 'All' not in allowed_strand_pair_tags:
            line = fp.readline()
            continue

        # Sanity check: Should never be entered, if analysis is restricted to AA interactions
        if tsss_a == '':
            print("\t\t[Warning] No TSS for first digest!")
            first_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # Sanity check: Should never be entered, if analysis is restricted to AA interactions
        if tsss_b == '':
            print("\t\t[Warning] No TSS for second digest!")
            second_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # Restrict analysis to specified set of categories for directed interactions (typicall we use both, i.e. 'S' and 'T')
        if interaction_category in allowed_interaction_categories_directed:

            # Increment conter for directed interactions
            directed_interaction_num += 1

            # Add digest regions to set for directed interactions
            directed_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            directed_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            # Add associated gene symbols to set for directed interactions
            for symbol in syms_a.split(","):
                directed_digests_symbols_set.add(symbol)
            for symbol in syms_b.split(","):
                directed_digests_symbols_set.add(symbol)

        # Restrict analysis to specified set of categories for undirected interactions (typicall we use 'URAA' only)
        if interaction_category in allowed_interaction_categories_undirected:

            # Increment conter for undirected interactions
            undirected_interaction_num += 1

            # Add digest regions to set for undirected interactions
            undirected_digests_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            undirected_digests_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            # Add associated gene symbols to set for directed interactions
            for symbol in syms_a.split(","):
                undirected_digests_symbols_set.add(symbol)
            for symbol in syms_b.split(","):
                undirected_digests_symbols_set.add(symbol)

        line = fp.readline()

# Sanity check: Report interactions with no TSS on one or both digests (should never happen, if analysis is restricted to AA interactions)
if 0 < first_digest_without_tss_num:
    print("\t[WARNING] There were " + str(first_digest_without_tss_num) + " interactions without promoters on the first digest!")
if 0 < second_digest_without_tss_num:
    print("\t[WARNING] There were " + str(second_digest_without_tss_num) + " interactions without promoters on the second digest!")

fp.close()

print("\t[INFO] ... done.")


### Get exclusive digests (D\U and U\D) and write coordinates for sequence analysis to BED files
################################################################################################

print("\t[INFO] Getting exclusive digests and writing coordinates for sequence analysis ...")

if comparative_sets == "D-VS-UWOD":
    directed_digests_comparative_set = directed_digests_set
elif comparative_sets == "DWOU-VS-UWOD":
    directed_digests_comparative_set = directed_digests_set.difference(undirected_digests_set)
else:
    print("\t\t[ERROR] Comparative sets should be either X or Y. Abort!")
    exit(1)

undirected_digests_comparative_set = undirected_digests_set.difference(directed_digests_set)
directed_intersect_undirected_digests_set = directed_digests_set.intersection(undirected_digests_set)
directed_union_directed_digests_set = directed_digests_set.union(undirected_digests_set)

cnt = 1 # use consecutive numbers as unique IDs
print("\t\t[INFO] Writing to file: " + bed_file_name_directed_digests)
for digest in directed_digests_comparative_set:
    bed_stream_name_directed_digests.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1
cnt = 1
print("\t\t[INFO] Writing to file: " + bed_file_name_undirected_digests)
for digest in undirected_digests_comparative_set:
    bed_stream_name_undirected_digests.write(digest + "\t" + str(cnt) + "\n")
    cnt += 1

bed_stream_name_directed_digests.close()
bed_stream_name_undirected_digests.close()

print("\t[INFO] ... done.")

print("[INFO] ... done with first pass.")


########################################################################################################################
########################################################################################################################

### Get unique TSS associated with unique exclusive digests (D\U and U\D)
#########################################################################

# Sets that will store coordinates and strands of TSS on digests
directed_digests_tss_comparative_set = set()
undirected_digests_tss_comparative_set = set()

# Arrays that will store interaction sizes
interaction_sizes_directed_digests_array = []
interaction_sizes_undirected_digests_array = []

print("[INFO] Second pass: Determining promoters on unique exclusive digests ...")

print("\t[INFO] Iterating interactions with gene symbols ...")

with gzip.open(interaction_gs_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline().rstrip()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse line from interaction file with gene symbols
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category = \
            parse_interaction_line_with_gene_symbols(line.rstrip())

        # Skip lines/interactions with the wrong enrichment pair tag (typically we use only AA)
        if enrichment_pair_tag not in allowed_enrichment_pair_tags: # use only interactions with specified digest enrichment pair tags (II, AI, AA)
            line = fp.readline()
            continue

        # Skrip line with wrong strand symbol pair tag (typically we use no restriction, i.e. 'All')
        if strand_pair_tag not in allowed_strand_pair_tags and 'All' not in allowed_strand_pair_tags: # use only interactions with specified digest strand pair tags (-/-, +/+, -/+, +/-, -/d, ...) or 'All'
            line = fp.readline()
            continue

        # Sanity check: Should never be entered, if analysis is restricted to AA interactions
        if tsss_a == '':
            print("\t\t[Warning] No TSS for first digest!")
            first_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # Sanity check: Should never be entered, if analysis is restricted to AA interactions
        if tsss_b == '':
            print("\t\t[Warning] No TSS for first digest!")
            second_digest_without_tss_num += 1
            line = fp.readline()
            continue

        # Get coordinates of the first digest from interaction line
        d1_coord = str(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))

        # Check if the first digest belongs to a directed interaction and does not additionally interact with a digest from the undirected reference interaction set
        if d1_coord in directed_digests_comparative_set:
            # Add all TSS positions on this digest to a set of TSS
            for tss_a in tsss_a.split(','): # ToDo: If we do want to avoid overlapping promoters, we could do it here
                directed_digests_tss_comparative_set.add(tss_a)
        # Check if the first digest belongs to a undirected reference interaction and does not additionally interact with a digest from the directed interaction set
        if d1_coord in undirected_digests_comparative_set:
            # Add all TSS positions on this digest to a set of TSS
            for tss_a in tsss_a.split(','): # ToDo: If we do want to avoid overlapping promoters, we could do it here
                undirected_digests_tss_comparative_set.add(tss_a)

        # Get coordinates of the second digest from interaction line
        d2_coord = str(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

        # Do the same as for the first digest
        if d2_coord in directed_digests_comparative_set:
            for tss_b in tsss_b.split(','): # ToDo: If we do want to avoid overlapping promoters, we could do it here
                directed_digests_tss_comparative_set.add(tss_b)
        if d2_coord in undirected_digests_comparative_set:
            for tss_b in tsss_b.split(','): # ToDo: If we do want to avoid overlapping promoters, we could do it here
                undirected_digests_tss_comparative_set.add(tss_b)

        # Get distance between directed interacting digests and strand pair tags, if both digests digest do not interact with a digests of the set of undirected reference interactions
        if interaction_category in allowed_interaction_categories_directed and d1_coord in directed_digests_comparative_set and d2_coord in directed_digests_comparative_set:

            # Collect distance between interacting digests
            interaction_sizes_directed_digests_array.append(sta_b - end_a)

            # Collect strand pair tag
            if strand_pair_tag in strand_pair_tag_directed_dict:
                strand_pair_tag_directed_dict[strand_pair_tag] += 1
            else:
                strand_pair_tag_directed_dict[strand_pair_tag] = 1

            # Count interaction partners within directed interactions for each digest
            if chr_a + "\t" + str(sta_a) + "\t" + str(end_a) in interaction_partners_per_digest_directed_dict.keys():
                interaction_partners_per_digest_directed_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] += 1
            else:
                interaction_partners_per_digest_directed_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] = 1

            if chr_b + "\t" + str(sta_b) + "\t" + str(end_b) in interaction_partners_per_digest_directed_dict.keys():
                interaction_partners_per_digest_directed_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] += 1
            else:
                interaction_partners_per_digest_directed_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] = 1

        # Get distance between undirected interacting digests and strand pair tags, if both digests digest do not interact with a digests of the set of directed interactions
        if interaction_category in allowed_interaction_categories_undirected and d1_coord in undirected_digests_comparative_set and d2_coord in undirected_digests_comparative_set:

            # Collect distance between interacting digests
            interaction_sizes_undirected_digests_array.append(sta_b - end_a)

            # Collect strand pair tag
            if strand_pair_tag in strand_pair_tag_undirected_dict:
                strand_pair_tag_undirected_dict[strand_pair_tag] += 1
            else:
                strand_pair_tag_undirected_dict[strand_pair_tag] = 1

            # Count interaction partners within undirected reference interactions (typically 'URAA') for each digest
            if chr_a + "\t" + str(sta_a) + "\t" + str(end_a) in interaction_partners_per_digest_undirected_dict.keys():
                interaction_partners_per_digest_undirected_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] += 1
            else:
                interaction_partners_per_digest_undirected_dict[chr_a + "\t" + str(sta_a) + "\t" + str(end_a)] = 1

            if chr_b + "\t" + str(sta_b) + "\t" + str(end_b) in interaction_partners_per_digest_undirected_dict.keys():
                interaction_partners_per_digest_undirected_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] += 1
            else:
                interaction_partners_per_digest_undirected_dict[chr_b + "\t" + str(sta_b) + "\t" + str(end_b)] = 1

        line = fp.readline()

print("\t[INFO] ... done.")


### Write promoter coordinates to BED files for sequence analysis
#################################################################

print("\t[INFO] Writing promoter coordinates to BED files for sequence analysis ...")

# Report names of files that are going to be created for TSS on directed digests
print("\t\t[INFO] Writing to file: " + bed_file_name_directed_digests_tss_minus)
print("\t\t[INFO] Writing to file: " + bed_file_name_directed_digests_tss_plus)

cnt_minus = 1 # use consecutive numbers as unique IDs
cnt_plus = 1
for tss in directed_digests_tss_comparative_set:

    # Split string with TSS coordinate and strand and add specified number of bases in up and downstream direction
    A = tss.split(':')
    chr = A[0]
    sta = int(A[1]) - up_dist
    if sta < 0: # avoid coordinates smaller than 0
        sta = 0
    end = int(A[1]) + down_dist
    if chrom_sizes_dict[chr] < end: # avoid coordinates greater than chromosome size
        end = chrom_sizes_dict[chr] - 1
    strand = A[2]

    # Write promoter coordinates to separate files for each strand
    if strand == '-':
        bed_stream_name_directed_digests_tss_minus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        bed_stream_name_directed_digests_tss_plus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
        cnt_plus += 1
    else:
        print("[Warning] Strand of TSS was neither \'-\' nor \'+\'!") # should never happen

# Close output streams
bed_stream_name_directed_digests_tss_minus.close()
bed_stream_name_directed_digests_tss_plus.close()
directed_digests_tss_comparative_set_size = len(directed_digests_tss_comparative_set)

# Report names of files that are going to be created for TSS on undirected digests
print("\t\t[INFO] Writing to file: " + bed_file_name_undirected_digests_tss_minus)
print("\t\t[INFO] Writing to file: " + bed_file_name_undirected_digests_tss_plus)

cnt_minus = 1
cnt_plus = 1
for tss in undirected_digests_tss_comparative_set:
    A = tss.split(':')
    chr = A[0]
    sta = int(A[1]) - up_dist
    if sta < 0:
        sta = 0
    end = int(A[1]) + down_dist
    if chrom_sizes_dict[chr] < end:
        end = chrom_sizes_dict[chr] - 1
    strand = A[2]
    if strand == '-':
        bed_stream_name_undirected_digests_tss_minus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_minus) + "\n")
        cnt_minus += 1
    elif strand == '+':
        bed_stream_name_undirected_digests_tss_plus.write(chr + "\t" + str(sta) + "\t" + str(end) + "\t" + str(cnt_plus) + "\n")
        cnt_plus += 1
    else:
        print("[Warning] Strand of TSS was neither \'-\' nor \'+\'!") # should never happen

# Close output streams
bed_stream_name_undirected_digests_tss_minus.close()
bed_stream_name_undirected_digests_tss_plus.close()
undirected_digests_tss_comparative_set_size = len(undirected_digests_tss_comparative_set)

print("\t[INFO] ... done.")

print("[INFO] ... done with second pass.")


########################################################################################################################
########################################################################################################################

### Calculate statistics on unique exclusive digests and associated promoters
#############################################################################

print("[INFO] Calculating statistics on unique exclusive digests and associated promoters ...")

### Determine and write digest sizes to file
############################################

print("\t[INFO] Determining and writing digest sizes ...")

tab_file_name_directed_digests_size_distribution = out_prefix + "_directed_digests_size_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_directed_digests_size_distribution)
directed_digest_q1, directed_digest_q2, directed_digest_q3, directed_digest_mean = determine_digest_sizes_and_write_to_file(tab_file_name_directed_digests_size_distribution, directed_digests_comparative_set)

tab_file_name_undirected_digests_size_distribution = out_prefix + "_undirected_digests_size_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_undirected_digests_size_distribution)
undirected_digest_q1, undirected_digest_q2, undirected_digest_q3, undirected_digest_mean = determine_digest_sizes_and_write_to_file(tab_file_name_undirected_digests_size_distribution, undirected_digests_comparative_set)
print("\t[INFO] ... done.")


### Write interaction sizes to file
###################################

print("\t[INFO] Writing interaction sizes to file ...")

tab_file_name_directed_interaction_size_distribution = out_prefix + "_directed_interaction_size_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_directed_interaction_size_distribution)
q1_interaction_sizes_directed_digests, q2_interaction_sizes_directed_digests, q3_interaction_sizes_directed_digests, mean_interaction_sizes_directed_digests = \
    write_interaction_sizes_to_file(interaction_sizes_directed_digests_array, tab_file_name_directed_interaction_size_distribution)

tab_file_name_undirected_interaction_size_distribution = out_prefix + "_undirected_interaction_size_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_undirected_interaction_size_distribution)
q1_interaction_sizes_undirected_digests, q2_interaction_sizes_undirected_digests, q3_interaction_sizes_undirected_digests, mean_interaction_sizes_undirected_digests = \
    write_interaction_sizes_to_file(interaction_sizes_undirected_digests_array, tab_file_name_undirected_interaction_size_distribution)

print("\t[INFO] ... done.")


### Write number of interaction partners per digest
###################################################

print("\t[INFO] Writing number of interaction partners per digest ...")

tab_file_name_directed_digests_interaction_partner_per_digest_distribution = out_prefix + "_directed_digests_interaction_partner_per_digest_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_directed_digests_interaction_partner_per_digest_distribution)
mean_number_of_interaction_partners_directed = write_number_of_interaction_partners_per_digest(tab_file_name_directed_digests_interaction_partner_per_digest_distribution, interaction_partners_per_digest_directed_dict)

tab_file_name_undirected_digests_interaction_partner_per_digest_distribution = out_prefix + "_undirected_digests_interaction_partner_per_digest_distribution.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_undirected_digests_interaction_partner_per_digest_distribution)
mean_number_of_interaction_partners_undirected = write_number_of_interaction_partners_per_digest(tab_file_name_undirected_digests_interaction_partner_per_digest_distribution, interaction_partners_per_digest_undirected_dict)

print("\t[INFO] ... done.")


### Get unique exclusive genes symbols for GO analysis
######################################################

print("\t[INFO] Getting unique exclusive genes symbols for GO analysis ...")

tab_file_name_directed_digest_symbols = out_prefix + "_directed_digest_symbols.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_directed_digest_symbols)

if comparative_sets == "D-VS-UWOD":
    empty_set = set() # do not subtract anything from directed digests
    symbols_abs_directed, symbols_znf_abs_directed, symbols_znf_rel_directed = \
        get_unique_exclusive_gene_symbols_and_write_to_file(directed_digests_symbols_set, empty_set, tab_file_name_directed_digest_symbols, "ZNF")
elif comparative_sets == "DWOU-VS-UWOD":
    symbols_abs_directed, symbols_znf_abs_directed, symbols_znf_rel_directed = \
        get_unique_exclusive_gene_symbols_and_write_to_file(directed_digests_symbols_set, undirected_digests_symbols_set, tab_file_name_directed_digest_symbols, "ZNF")
else:
    print("\t\t[ERROR] Comparative sets should be either D-VS-UWOD or DWOU-VS-UWOD. Abort!")
    exit(1)

tab_file_name_undirected_digest_symbols = out_prefix + "_undirected_digest_symbols.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_undirected_digest_symbols)
symbols_abs_undirected, symbols_znf_abs_undirected, symbols_znf_rel_undirected = \
    get_unique_exclusive_gene_symbols_and_write_to_file(undirected_digests_symbols_set, directed_digests_symbols_set, tab_file_name_undirected_digest_symbols, "ZNF")

print("\t[INFO] ... done.")


### Convert BED with unique exclusive digests to FASTA
######################################################

print("\t[INFO] Converting BED with unique exclusive digests to FASTA ...")

# get FASTA files with DIGEST sequences of directed interactions for calculation of sequence statistics
fasta_file_name_directed_digests = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_directed_digests)

# get FASTA files with DIGEST sequences of undirected interactions for calculation of sequence statistics
fasta_file_name_undirected_digests = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_undirected_digests)


# get FASTA files with sequences around promoters on '-' strand that are associated with directed interactions for CentriMo analysis
fasta_file_name_directed_digests_tss_minus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_minus)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_directed_digests_tss_minus)

# get FASTA files with sequences around promoters on '+' strand that are associated with directed interactions for CentriMo analysis
fasta_file_name_directed_digests_tss_plus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_plus)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_directed_digests_tss_plus)

# get FASTA files with sequences around promoters on '-' strand that are associated with undirected interactions for CentriMo analysis
fasta_file_name_undirected_digests_tss_minus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_minus)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_undirected_digests_tss_minus)

# get FASTA files with sequences around promoters on '+' strand that are associated with undirected interactions for CentriMo analysis
fasta_file_name_undirected_digests_tss_plus = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_plus)
print("\t\t[INFO] Wrote to file: " + fasta_file_name_undirected_digests_tss_plus)

# concat BED files with promoter regions for plus and minus strand
sys_cmd = 'cat ' + bed_file_name_directed_digests_tss_minus + ' ' + bed_file_name_directed_digests_tss_plus + ' > ' + bed_file_name_directed_digests_tss
os.system(sys_cmd)
sys_cmd = 'cat ' + bed_file_name_undirected_digests_tss_minus + ' ' + bed_file_name_undirected_digests_tss_plus + ' > ' + bed_file_name_undirected_digests_tss
os.system(sys_cmd)

# remove separate BED files for different strands
sys_cmd = 'rm ' + bed_file_name_directed_digests_tss_minus + ' ' + bed_file_name_directed_digests_tss_plus
os.system(sys_cmd)
sys_cmd = 'rm ' + bed_file_name_undirected_digests_tss_minus + ' ' + bed_file_name_undirected_digests_tss_plus
os.system(sys_cmd)

# merge promoter regions
bed_file_name_directed_digests_tss_merged = merge_bed(bedtools_path, bed_file_name_directed_digests_tss)
print("\t\t[INFO] Wrote to file: " + bed_file_name_directed_digests_tss_merged)
bed_file_name_undirected_digests_tss_merged = merge_bed(bedtools_path, bed_file_name_undirected_digests_tss)
print("\t\t[INFO] Wrote to file: " + bed_file_name_undirected_digests_tss_merged)

# get FASTA with merged promoter sequences for calculation of sequence statistics
fasta_file_name_directed_digests_tss_merged = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_merged) # used for sequence statistics of promoters
print("\t\t[INFO] Wrote to file: " + fasta_file_name_directed_digests_tss_merged)
fasta_file_name_undirected_digests_tss_merged = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_merged) # used for sequence statistics of promoters
print("\t\t[INFO] Wrote to file: " + fasta_file_name_undirected_digests_tss_merged)

# trim merged promoter regions to a uniform length for motif analysis with DREME
target_region_size_center_trimmed_directed, total_region_num_center_trimmed_directed, trimmed_region_num_frac_center_trimmed_directed, avg_region_size_center_trimmed_directed, region_max_size_center_trimmed_directed, total_trimmed_region_size_frac_center_trimmed_directed, bed_file_name_directed_digests_tss_merged_center_trimmed =\
    center_trim_bed(bed_file_name_directed_digests_tss_merged, up_dist, down_dist)
target_region_size_center_trimmed_undirected, total_region_num_center_trimmed_undirected, trimmed_region_num_frac_center_trimmed_undirected, avg_region_size_center_trimmed_undirected, region_max_size_center_trimmed_undirected, total_trimmed_region_size_frac_center_trimmed_undirected, bed_file_name_undirected_digests_tss_merged_center_trimmed =\
    center_trim_bed(bed_file_name_undirected_digests_tss_merged, up_dist, down_dist)

# convert merged and trimmed promoter regions to FASTA
fasta_file_name_directed_digests_tss_merged_center_trimmed = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_directed_digests_tss_merged_center_trimmed) # could be removed after repeat masking
fasta_file_name_directed_digests_tss_merged_center_trimmed_masked = mask_repeats(fasta_file_name_directed_digests_tss_merged_center_trimmed)
fasta_file_name_undirected_digests_tss_merged_center_trimmed = convert_bed_to_fasta(bedtools_path, genome_fasta_path, bed_file_name_undirected_digests_tss_merged_center_trimmed) # could be removed after repeat masking
fasta_file_name_undirected_digests_tss_merged_center_trimmed_masked = mask_repeats(fasta_file_name_undirected_digests_tss_merged_center_trimmed)

# remove unmasked FASTA files with sequences of merged, centered and trimmed promoters
sys_cmd = 'rm ' + fasta_file_name_directed_digests_tss_merged_center_trimmed + ' ' + fasta_file_name_undirected_digests_tss_merged_center_trimmed
os.system(sys_cmd)

print("\t[INFO] ... done.")


### Determine base frequencies and write to file
################################################

print("\t[INFO] Determining base frequencies and writing to file ...")

tab_file_name_base_frequencies = out_prefix + "_base_frequencies.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_base_frequencies)
tab_stream_name_base_frequencies = open(tab_file_name_base_frequencies, 'wt')

# Directed digests
header_line, value_line, repeat_content_directed_digests, gc_content_repeat_directed_digests, gc_content_non_repeat_directed_digests, gc_content_total_directed_digests = \
    get_base_frequencies(fasta_file_name_directed_digests)
tab_stream_name_base_frequencies.write(header_line + "\n")
tab_stream_name_base_frequencies.write(value_line + "\n")

# Undirected digests
header_line, value_line, repeat_content_undirected_digests, gc_content_repeat_undirected_digests, gc_content_non_repeat_undirected_digests, gc_content_total_undirected_digests = \
    get_base_frequencies(fasta_file_name_undirected_digests)
tab_stream_name_base_frequencies.write(value_line + "\n")

# Merged promoters on directed digests
header_line, value_line, repeat_content_directed_digests_tss_merged, gc_content_repeat_directed_digests_tss_merged, gc_content_non_repeat_directed_digests_tss_merged, gc_content_total_directed_digests_tss_merged = \
    get_base_frequencies(fasta_file_name_directed_digests_tss_merged)
tab_stream_name_base_frequencies.write(value_line + "\n")

# Merged promoters on undirected digests
header_line, value_line, repeat_content_undirected_digests_tss_merged, gc_content_repeat_undirected_digests_tss_merged, gc_content_non_repeat_undirected_digests_tss_merged, gc_content_total_undirected_digests_tss_merged = \
    get_base_frequencies(fasta_file_name_undirected_digests_tss_merged)
tab_stream_name_base_frequencies.write(value_line + "\n")

tab_stream_name_base_frequencies.close()

# remove FASTA files for calculation of sequence statistics on digests and promoters
sys_cmd = 'rm ' + fasta_file_name_directed_digests + ' ' + fasta_file_name_undirected_digests
os.system(sys_cmd)
sys_cmd = 'rm ' + fasta_file_name_directed_digests_tss_merged + ' ' + fasta_file_name_undirected_digests_tss_merged
os.system(sys_cmd)

print("\t[INFO] ... done.")

### Determine base frequencies and write to file
################################################

print("\t[INFO] Determining frequencies of strand pair tags and writing to file ...")

tab_file_name_strand_pair_tag_frequencies = out_prefix + "_strand_pair_tag_frequencies.tab"
print("\t\t[INFO] Writing to file: " + tab_file_name_strand_pair_tag_frequencies)
tab_stream_name_strand_pair_tag_frequencies = open(tab_file_name_strand_pair_tag_frequencies, 'wt')

# Get total numbers of strand pair tags
strand_pair_tag_directed_num = 0
strand_pair_tag_undirected_num = 0
for key in strand_pair_tag_list:
    strand_pair_tag_directed_num = strand_pair_tag_directed_num + strand_pair_tag_directed_dict[key]
    strand_pair_tag_undirected_num = strand_pair_tag_undirected_num + strand_pair_tag_undirected_dict[key]

# Write frequencies to file
for key in strand_pair_tag_list:
    tab_stream_name_strand_pair_tag_frequencies.write('\'' + key + '\'' + "\t" +
                                                      "{:.2f}".format(100 * strand_pair_tag_directed_dict[key] / strand_pair_tag_directed_num) + "\t" +
                                                      "{:.2f}".format(100 * strand_pair_tag_undirected_dict[key] / strand_pair_tag_undirected_num) + "\t" +
                                                      str(strand_pair_tag_directed_dict[key]) + "\t" +
                                                      str(strand_pair_tag_undirected_dict[key]) +
                                                      "\n")
tab_stream_name_strand_pair_tag_frequencies.close()

print("\t[INFO] ... done.")

print("[INFO] ... done with calculation of statistics on unique exclusive digests and associated promoters.")


### Create file with interaction and digest statistics
######################################################

print("[INFO] Creating file with interaction and digest statistics ...")
print("\t[INFO] Number of directed interactions: " + str(directed_interaction_num))
print("\t[INFO] Number of undirected interactions: " + str(undirected_interaction_num))

print("\t[INFO] Number of all unique directed digests: " + str(len(directed_digests_set)))
print("\t[INFO] Number of all unique undirected digests: " + str(len(undirected_digests_set)))

print("\t[INFO] Number of unique directed digests that are not undirected: " + str(len(directed_digests_comparative_set)))
print("\t[INFO] Number of unique undirected digests that are not directed: " + str(len(undirected_digests_comparative_set)))

print("\t[INFO] Number of intersecting digests: " + str(len(directed_intersect_undirected_digests_set)))
print("\t[INFO] Size of the union of directed and undirected digests: " + str(len(directed_union_directed_digests_set)))

connectivity_factor_directed = 1-len(directed_digests_set)/(directed_interaction_num*2) # zero, if no digest belongs to more than one interaction
print("\t[INFO] Connectivity factor of directed digests: " + str("{:.2f}".format(connectivity_factor_directed)))

connectivity_factor_undirected = 1-len(undirected_digests_set)/(undirected_interaction_num*2)
print("\t[INFO] Connectivity factor of undirected digests: " + str("{:.2f}".format(connectivity_factor_undirected)))

# Jaccard index
jaccard_index = len(directed_intersect_undirected_digests_set) / len(directed_union_directed_digests_set)
print("\t[INFO] Jaccard index: " + str("{:.2f}".format(jaccard_index)))

# Szymkiewicz–Simpson coefficient
szymkiewicz_simpson_coefficient = len(directed_intersect_undirected_digests_set) / (min(len(directed_digests_set), len(undirected_digests_set)))
print("\t[INFO] Szymkiewicz–Simpson coefficient: " + str("{:.2f}".format(szymkiewicz_simpson_coefficient)))

print("\t[INFO] Repeat content of directed digests: " + str(repeat_content_directed_digests))
print("\t[INFO] Repeat content of undirected digests: " + str(repeat_content_undirected_digests))
print("\t[INFO] Repeat content of promoters on directed digests (merged): " + str(repeat_content_directed_digests_tss_merged))
print("\t[INFO] Repeat content of promoters on undirected digests (merged): " + str(repeat_content_undirected_digests_tss_merged))

print("\t[INFO] GC content within repeat regions of directed digests: " + str(gc_content_repeat_directed_digests))
print("\t[INFO] GC content within repeat regions of undirected digests: " + str(gc_content_repeat_undirected_digests))
print("\t[INFO] GC content within repeat regions of promoters on directed digests (merged): " + str(gc_content_repeat_directed_digests_tss_merged))
print("\t[INFO] GC content within repeat regions of promoters on undirected digests (merged): " + str(gc_content_repeat_undirected_digests_tss_merged))

print("\t[INFO] GC content within non repeat regions of directed digests: " + str(gc_content_non_repeat_directed_digests))
print("\t[INFO] GC content within non repeat regions of undirected digests: " + str(gc_content_non_repeat_undirected_digests))
print("\t[INFO] GC content within non repeat regions of promoters on directed digests (merged): " + str(gc_content_non_repeat_directed_digests_tss_merged))
print("\t[INFO] GC content within non repeat regions of promoters on undirected digests (merged): " + str(gc_content_non_repeat_undirected_digests_tss_merged))

print("\t[INFO] Total GC content of directed digests: " + str(gc_content_total_directed_digests))
print("\t[INFO] Total GC content of undirected digests: " + str(gc_content_total_undirected_digests))
print("\t[INFO] Total GC content of promoters on directed digests (merged): " + str(gc_content_total_directed_digests_tss_merged))
print("\t[INFO] Total GC content of promoters on undirected digests (merged): " + str(gc_content_total_undirected_digests_tss_merged))

print("\t[INFO] Size distribution of directed digests")
print("\t\t[INFO] Directed digests")
print("\t\t\t[INFO] Q1: " + str(directed_digest_q1) + " bp")
print("\t\t\t[INFO] Q2: " + str(directed_digest_q2) + " bp (Median)")
print("\t\t\t[INFO] Q3: " + str(directed_digest_q3) + " bp")
print("\t\t\t[INFO] Mean: " + str(directed_digest_mean) + " bp")
print("\t\t[INFO] Undirected digests")
print("\t\t\t[INFO] Q1: " + str(undirected_digest_q1) + " bp")
print("\t\t\t[INFO] Q2: " + str(undirected_digest_q2) + " bp (Median)")
print("\t\t\t[INFO] Q3: " + str(undirected_digest_q3) + " bp")
print("\t\t\t[INFO] Mean: " + str(undirected_digest_mean) + " bp")
print("\t[INFO] Mean number of interaction partners of directed digests: " + mean_number_of_interaction_partners_directed)
print("\t[INFO] Mean number of interaction partners of undirected digests: " + mean_number_of_interaction_partners_undirected)

print("\t[INFO] Unique exclusive gene symbols")
print("\t\t[INFO] Directed digests")
print("\t\t\t[INFO] Total number exclusive unique gene symbols for directed interactions: " + str(symbols_abs_directed))
print("\t\t\t[INFO] Number exclusive unique gene symbols with \'ZNF\' substring for directed interactions: " + str(symbols_znf_abs_directed))
print("\t\t\t[INFO] Fraction exclusive unique gene symbols with \'ZNF\' substring for directed interactions: " + str(symbols_znf_rel_directed))
print("\t\t[INFO] Undirected digests")
print("\t\t\t[INFO] Total number exclusive unique gene symbols for undirected interactions: " + str(symbols_abs_undirected))
print("\t\t\t[INFO] Number exclusive unique gene symbols with \'ZNF\' substring for undirected interactions: " + str(symbols_znf_abs_undirected))
print("\t\t\t[INFO] Fraction exclusive unique gene symbols with \'ZNF\' substring for undirected interactions: " + str(symbols_znf_rel_undirected))
print("\t[INFO] Interaction sizes")
print("\t\t[INFO] Directed digests")
print("\t\t\t[INFO] Q1: " + str(q1_interaction_sizes_directed_digests))
print("\t\t\t[INFO] Q2: " + str(q2_interaction_sizes_directed_digests))
print("\t\t\t[INFO] Q3: " + str(q3_interaction_sizes_directed_digests))
print("\t\t\t[INFO] Mean: " + str(mean_interaction_sizes_directed_digests))
print("\t\t[INFO] Undirected digests")
print("\t\t\t[INFO] Q1: " + str(q1_interaction_sizes_undirected_digests))
print("\t\t\t[INFO] Q2: " + str(q2_interaction_sizes_undirected_digests))
print("\t\t\t[INFO] Q3: " + str(q3_interaction_sizes_undirected_digests))
print("\t\t\t[INFO] Mean: " + str(mean_interaction_sizes_undirected_digests))

print("\t[INFO] Strand pair tags")
print("\t\tstrand_pair_tag\tfrac_directed\tfrac_undirected\tabs_directed\tabs_undirected")
for key in strand_pair_tag_list:
    print("\t\t" + '\'' + key + '\'' + "\t" +
          "{:.2f}".format(100 * strand_pair_tag_directed_dict[key] / strand_pair_tag_directed_num) + "\t" + "{:.2f}".format(100 * strand_pair_tag_undirected_dict[key] / strand_pair_tag_undirected_num) + "\t" +
          str(strand_pair_tag_directed_dict[key]) + "\t" + str(strand_pair_tag_undirected_dict[key])
          )

print("\t[INFO] Merged and trimmed promoters")
print("\t\t[INFO] Directed")
print("\t\t\t[INFO] Target region size: " + str(target_region_size_center_trimmed_directed))
print("\t\t\t[INFO] Total number of input regions (merged): " + str(total_region_num_center_trimmed_directed))
print("\t\t\t[INFO] Fraction of trimmed regions: " + str(trimmed_region_num_frac_center_trimmed_directed) + "%")
print("\t\t\t[INFO] Average size of regions: " + str(avg_region_size_center_trimmed_directed))
print("\t\t\t[INFO] Size of largest region: " + str(region_max_size_center_trimmed_directed))
print("\t\t\t[INFO] Fraction of trimmed sequence: " + str(total_trimmed_region_size_frac_center_trimmed_directed) + "%")
print("\t\t\t[INFO] BED file with trimmed regions: " + bed_file_name_directed_digests_tss_merged_center_trimmed)

print("\t\t[INFO] Unirected")
print("\t\t\t[INFO] Target region size: " + str(target_region_size_center_trimmed_undirected))
print("\t\t\t[INFO] Total number of input regions (merged): " + str(total_region_num_center_trimmed_undirected))
print("\t\t\t[INFO] Fraction of trimmed regions: " + str(trimmed_region_num_frac_center_trimmed_undirected) + "%")
print("\t\t\t[INFO] Average size of regions: " + str(avg_region_size_center_trimmed_undirected))
print("\t\t\t[INFO] Size of largest region: " + str(region_max_size_center_trimmed_undirected))
print("\t\t\t[INFO] Fraction of trimmed sequence: " + str(total_trimmed_region_size_frac_center_trimmed_undirected) + "%")
print("\t\t\t[INFO] BED file with trimmed regions: " + bed_file_name_undirected_digests_tss_merged_center_trimmed)


tab_file_name_interaction_and_digest_statistics = out_prefix + "_interaction_and_digest_statistics.tab"
tab_file_stream_interaction_and_digest_statistics = open(tab_file_name_interaction_and_digest_statistics, 'wt')

tab_file_stream_interaction_and_digest_statistics.write(
    "out_prefix" + "\t" +

    "directed_interaction_num" + "\t" +
    "undirected_interaction_num" + "\t" +

    "directed_digests_set_size" + "\t" +
    "undirected_digests_set_size" + "\t" +

    "directed_digests_comparative_set_size" + "\t" +
    "undirected_digests_comparative_set_size" + "\t" +

    "directed_intersect_undirected_digests_set" + "\t" +
    "directed_union_directed_digests_set_size" + "\t" +

    "connectivity_factor_directed" + "\t" +
    "connectivity_factor_undirected" + "\t" +

    "jaccard_index" + "\t" +
    "szymkiewicz_simpson_coefficient" + "\t" +

    "repeat_content_directed_digests" + "\t" +
    "repeat_content_undirected_digests" + "\t" +
    "repeat_content_directed_digests_tss_merged" + "\t" +
    "repeat_content_undirected_digests_tss_merged" + "\t" +

    "gc_content_repeat_directed_digests" + "\t" +
    "gc_content_repeat_undirected_digests" + "\t" +
    "gc_content_repeat_directed_digests_tss_merged" + "\t" +
    "gc_content_repeat_undirected_digests_tss_merged" + "\t" +

    "gc_content_non_repeat_directed_digests" + "\t" +
    "gc_content_non_repeat_undirected_digests" + "\t" +
    "gc_content_non_repeat_directed_digests_tss_merged" + "\t" +
    "gc_content_non_repeat_undirected_digests_tss_merged" + "\t" +

    "gc_content_total_directed_digests" + "\t" +
    "gc_content_total_undirected_digests" + "\t" +
    "gc_content_total_directed_digests_tss_merged" + "\t" +
    "gc_content_non_repeat_undirected_digests_tss_merged" + "\t" +

    "directed_digest_q1" + "\t" +
    "directed_digest_q2" + "\t" +
    "directed_digest_q3" + "\t" +
    "directed_digest_mean" + "\t" +

    "undirected_digest_q1" + "\t" +
    "undirected_digest_q2" + "\t" +
    "undirected_digest_q3" + "\t" +
    "undirected_digest_mean" + "\t" +

    "mean_number_of_interaction_partners_directed" + "\t" +
    "mean_number_of_interaction_partners_undirected" + "\t" +

    "symbols_abs_directed" + "\t" +
    "symbols_znf_abs_directed" + "\t" +
    "symbols_znf_rel_directed" + "\t" +

    "symbols_abs_undirected" + "\t" +
    "symbols_znf_abs_undirected" + "\t" +
    "symbols_znf_rel_undirected" + "\t" +

    "q1_interaction_sizes_directed_digests" + "\t" +
    "q2_interaction_sizes_directed_digests" + "\t" +
    "q3_interaction_sizes_directed_digests" + "\t" +
    "mean_interaction_sizes_directed_digests" + "\t" +

    "q1_interaction_sizes_undirected_digests" + "\t" +
    "q2_interaction_sizes_undirected_digests" + "\t" +
    "q3_interaction_sizes_undirected_digests" + "\t" +
    "mean_interaction_sizes_undirected_digests" + "\t" +

    "directed_digests_tss_comparative_set_size" + "\t" + # redundant promoters

    "total_region_num_center_trimmed_directed" + "\t" +  # merged promoters
    "target_region_size_center_trimmed_directed" + "\t" +
    "trimmed_region_num_frac_center_trimmed_directed" + "\t" +
    "avg_region_size_center_trimmed_directed" + "\t" +
    "total_trimmed_region_size_frac_center_trimmed_directed" + "\t" +

    "undirected_digests_tss_comparative_set_size" + "\t" + # redundant promoters

    "total_region_num_center_trimmed_undirected" + "\t" +  # merged promoters
    "target_region_size_center_trimmed_undirected" + "\t" +
    "trimmed_region_num_frac_center_trimmed_undirected" + "\t" +
    "avg_region_size_center_trimmed_undirected" + "\t" +
    "total_trimmed_region_size_frac_center_trimmed_undirected" +

    "\n")


# Write line with corresponding values
tab_file_stream_interaction_and_digest_statistics.write(

    out_prefix + "\t" +

    str(directed_interaction_num) + "\t" +
    str(undirected_interaction_num) + "\t" +

    str(len(directed_digests_set)) + "\t" +
    str(len(undirected_digests_set)) + "\t" +

    str(len(directed_digests_comparative_set)) + "\t" +
    str(len(undirected_digests_comparative_set)) + "\t" +

    str(len(directed_intersect_undirected_digests_set)) + "\t" +
    str(len(directed_union_directed_digests_set)) + "\t" +

    str("{:.2f}".format(connectivity_factor_directed)) + "\t" +
    str("{:.2f}".format(connectivity_factor_undirected)) + "\t" +

    str("{:.2f}".format(jaccard_index)) + "\t" +
    str("{:.2f}".format(szymkiewicz_simpson_coefficient)) + "\t" +

    str(repeat_content_directed_digests) + "\t" +
    str(repeat_content_undirected_digests) + "\t" +

    str(repeat_content_directed_digests_tss_merged) + "\t" +
    str(repeat_content_undirected_digests_tss_merged) + "\t" +

    str(gc_content_repeat_directed_digests) + "\t" +
    str(gc_content_repeat_undirected_digests) + "\t" +

    str(gc_content_repeat_directed_digests_tss_merged) + "\t" +
    str(gc_content_repeat_undirected_digests_tss_merged) + "\t" +

    str(gc_content_non_repeat_directed_digests) + "\t" +
    str(gc_content_non_repeat_undirected_digests) + "\t" +

    str(gc_content_non_repeat_directed_digests_tss_merged) + "\t" +
    str(gc_content_non_repeat_undirected_digests_tss_merged) + "\t" +

    str(gc_content_total_directed_digests) + "\t" +
    str(gc_content_total_undirected_digests) + "\t" +

    str(gc_content_total_directed_digests_tss_merged) + "\t" +
    str(gc_content_total_undirected_digests_tss_merged) + "\t" +

    str(directed_digest_q1) + "\t" +
    str(directed_digest_q2) + "\t" +
    str(directed_digest_q3) + "\t" +
    str(directed_digest_mean) + "\t" +

    str(undirected_digest_q1) + "\t" +
    str(undirected_digest_q2) + "\t" +
    str(undirected_digest_q3) + "\t" +
    str(undirected_digest_mean) + "\t" +

    str(mean_number_of_interaction_partners_directed) + "\t" +
    str(mean_number_of_interaction_partners_undirected) + "\t" +

    str(symbols_abs_directed) + "\t" +
    str(symbols_znf_abs_directed) + "\t" +
    str(symbols_znf_rel_directed) + "\t" +

    str(symbols_abs_undirected) + "\t" +
    str(symbols_znf_abs_undirected) + "\t" +
    str(symbols_znf_rel_undirected) + "\t" +

    str(q1_interaction_sizes_directed_digests) + "\t" +
    str(q2_interaction_sizes_directed_digests) + "\t" +
    str(q3_interaction_sizes_directed_digests) + "\t" +
    str(mean_interaction_sizes_directed_digests) + "\t" +

    str(q1_interaction_sizes_undirected_digests) + "\t" +
    str(q2_interaction_sizes_undirected_digests) + "\t" +
    str(q3_interaction_sizes_undirected_digests) + "\t" +
    str(mean_interaction_sizes_undirected_digests) + "\t" +

    str(directed_digests_tss_comparative_set_size) + "\t" +  # redundant promoters

    str(total_region_num_center_trimmed_directed) + "\t" +  # merged promoters
    str(target_region_size_center_trimmed_directed) + "\t" +
    str(trimmed_region_num_frac_center_trimmed_directed) + "\t" +
    str(avg_region_size_center_trimmed_directed) + "\t" +
    str(total_trimmed_region_size_frac_center_trimmed_directed) + "\t" +

    str(undirected_digests_tss_comparative_set_size) + "\t" +  # redundant promoters

    str(total_region_num_center_trimmed_undirected) + "\t" +  # merged promoters
    str(target_region_size_center_trimmed_undirected) + "\t" +
    str(trimmed_region_num_frac_center_trimmed_undirected) + "\t" +
    str(avg_region_size_center_trimmed_undirected) + "\t" +
    str(total_trimmed_region_size_frac_center_trimmed_undirected) +

    "\n")

tab_file_stream_interaction_and_digest_statistics.close()

print("[INFO] ... done.")
