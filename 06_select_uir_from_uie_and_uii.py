#!/usr/bin/env python

"""
This script takes an enhanced interaction file (created with the script '05_define_di_uie_and_uii.py') with the
following interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions
      i. UIE - Exclusive undirected interactions. Subset of UI. No digest interacts with a digest involved in DI
      ii. UII - Inclusive undirected interactions. Subset of UI. UI without UIE

that are specified in the third column.

Undirected reference interactions (UIR) that are comparable to DI with respect to enrichment status of digests, read
pair numbers and, optionally, interaction distances are selected from UIE and UII.

The script implements a two pass approach:

   1. Iterate interactions and collect information about DI with respect to enrichment status of digests, read pair
   numbers and interaction distances.
   2. Iterate interactions a second time and use the collected information to select reference interactions that are
   comparable to DI.

The script implements two different approaches for the selection of undirected reference interactions:

   1. Determine Q1 and Q3 for the distribution of read pair numbers within DI during the first pass and select all
   undirected interactions with read pair numbers between Q1 and Q2 during the second pass.
   2. Determine the exact distribution of read pair numbers within DI during the first pass and select the same number
   of undirected interactions for each read number during the second pass.

Use '--selection-approach q13' for the first and '--selection-approach exact' for the second approach.

Optionally, the selected reference interaction can be adjusted to DI with respect to interaction size (in addition to
the adjustment to read pair numbers). For this purpose, the distribution of distances within DI is determined during the
first pass and only undirected interactions with a size between Q1 and Q3 are selected as reference interactions. Use
'--adjust-to-interaction-distance' to adjust reference interaction to interaction size.

Directed interactions and undirected reference interactions are written to a file in enhanced interaction format with
the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_uir.tsv.gz'

Column 3 of the created file contains the following interaction categories:

   1. DIAA - Directed interactions within AA
   2. DIAI - Directed interactions within AI
   3. DIII - Directed interactions within II
   4. UIRAA - Undirected reference interactions within AA
   5. UIRAI - Undirected reference interactions within AI
   6. UIRII - Undirected reference interactions within II

Finally, summary statistics about numbers of interactions within different enrichment state categories (AA, AI and II)
and corresponding read pair numbers are printed to the screen and two boxplots are created, one for read pair numbers
and another one for interaction distances.

   '<OUT_PREFIX>_read_pair_number_boxplot.pdf'
   '<OUT_PREFIX>_interaction_distance_boxplot.pdf'
"""


import argparse
import gzip
import diachrscripts_toolkit
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


### Parse command line
######################

parser = argparse.ArgumentParser(description='Select reference interactions from exclusive undirected interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file created with \'05_define_di_uie_and_uii.py\'. Column 3 contains either \'DI\' or\'UIE\'.', required=True)
parser.add_argument('--selection-approach', help='Use \'q13\' to select undirected reference interactions with read pair numbers between Q1 and Q3 of the distribution for DI or \'exact\' to select a subset of undirected interactions with the almost same distribution of read pairs as compared to DI.', default='exact')
parser.add_argument('--adjust-to-interaction-distance', help='Adjust selected reference interactions to interaction distances of directed interactions.', action='store_true', default=False)
args = parser.parse_args()

out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
selection_approach = args.selection_approach
adjust_to_interaction_distance = args.adjust_to_interaction_distance

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --selection-approach: " + selection_approach)
print("\t[INFO] --adjust-to-interaction-distance: " + str(adjust_to_interaction_distance))

if not(selection_approach == 'q13' or selection_approach == 'exact'):
    print("\t[ERROR] --selection-approach must be \'q13\' or \'exact\'!")
    exit(1)


### Prepare variables, data structures and streams for output files
###################################################################

# Directed interactions
#######################

# Total number of interactions
dir_inter_num = 0

# Numbers of interactions within AA, AI and II
dir_inter_aa_num = 0
dir_inter_ai_num = 0
dir_inter_ii_num = 0

# Percentages of interactions within AA, AI and II
dir_inter_aa_percentage = None
dir_inter_ai_percentage = None
dir_inter_ii_percentage = None

dir_inter_aa_rp_q1 = None
dir_inter_aa_rp_median = None
dir_inter_aa_rp_q3 = None

# Arrays for read pair numbers within AA, AI and II
dir_inter_aa_rp_array = []
dir_inter_ai_rp_array = []
dir_inter_ii_rp_array = []

# First quantiles of read pair number distribution within AA, AI and II
dir_inter_aa_rp_q1 = None
dir_inter_ai_rp_q1 = None
dir_inter_ii_rp_q1 = None

# Second quantiles of read pair number distribution within AA, AI and II
dir_inter_aa_rp_median = None
dir_inter_ai_rp_median = None
dir_inter_ii_rp_median = None

# Third quantiles of read pair number distribution within AA, AI and II
dir_inter_aa_rp_q3 = None
dir_inter_ai_rp_q3 = None
dir_inter_ii_rp_q3 = None

# Arrays for interaction distances within AA, AI and II
dir_inter_aa_dist_array = []
dir_inter_ai_dist_array = []
dir_inter_ii_dist_array = []

# Undirected interactions
#########################

# Total number of interactions
undir_inter_num = 0

# Numbers of interactions within AA, AI and II
undir_inter_aa_num = 0
undir_inter_ai_num = 0
undir_inter_ii_num = 0

# Percentages of interactions within AA, AI and II
undir_inter_aa_percentage = None
undir_inter_ai_percentage = None
undir_inter_ii_percentage = None

# Arrays for read pair numbers within AA, AI and II
undir_inter_aa_rp_array = []
undir_inter_ai_rp_array = []
undir_inter_ii_rp_array = []

# Second quantiles of read pair number distribution within AA, AI and II
undir_inter_aa_rp_median = None
undir_inter_ai_rp_median = None
undir_inter_ii_rp_median = None

# Arrays for interaction distances within AA, AI and II
undir_inter_aa_dist_array = []
undir_inter_ai_dist_array = []
undir_inter_ii_dist_array = []

# Undirected reference interactions 1 ('q13' selection approach)
################################################################

# Total number of interactions
undir_ref_1_inter_num = 0

# Numbers of interactions within AA, AI and II
undir_ref_1_inter_aa_num = 0
undir_ref_1_inter_ai_num = 0
undir_ref_1_inter_ii_num = 0

# Percentages of interactions within AA, AI and II
undir_ref_1_inter_aa_percentage = None
undir_ref_1_inter_ai_percentage = None
undir_ref_1_inter_ii_percentage = None

# Arrays for read pair numbers within AA, AI and II
undir_ref_1_inter_aa_rp_array = []
undir_ref_1_inter_ai_rp_array = []
undir_ref_1_inter_ii_rp_array = []

# Second quantiles of read pair number distribution within AA, AI and II
undir_ref_1_inter_aa_rp_median = None
undir_ref_1_inter_ai_rp_median = None
undir_ref_1_inter_ii_rp_median = None

# Arrays for interaction distances within AA, AI and II
undir_ref_1_inter_aa_dist_array = []
undir_ref_1_inter_ai_dist_array = []
undir_ref_1_inter_ii_dist_array = []

# Undirected reference interactions 2 ('exact' selection approach)
##################################################################

# Total number of interactions
undir_ref_2_inter_num = 0

# Numbers of interactions within AA, AI and II
undir_ref_2_inter_aa_num = 0
undir_ref_2_inter_ai_num = 0
undir_ref_2_inter_ii_num = 0

# Percentages of interactions within AA, AI and II
undir_ref_2_inter_aa_percentage = None
undir_ref_2_inter_ai_percentage = None
undir_ref_2_inter_ii_percentage = None

# Arrays for read pair numbers within AA, AI and II
undir_ref_2_inter_aa_rp_array = []
undir_ref_2_inter_ai_rp_array = []
undir_ref_2_inter_ii_rp_array = []

# Second quantiles of read pair number distribution within AA, AI and II
undir_ref_2_inter_aa_rp_median = None
undir_ref_2_inter_ai_rp_median = None
undir_ref_2_inter_ii_rp_median = None

# Arrays for interaction distances within AA, AI and II
undir_ref_2_inter_aa_dist_array = []
undir_ref_2_inter_ai_dist_array = []
undir_ref_2_inter_ii_dist_array = []

# Dictionaries for read pair numbers of directed interactions within AA, AI and II
rp_dict_aa = {}
rp_dict_ai = {}
rp_dict_ii = {}

# Output files
##############

# Enhanced interaction file with interaction categories: DI, UIRAA, UIRAI and UIRII
enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_with_di_and_uir.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')

# PDF file with boxplots for distributions of read pair numbers
pdf_name_boxplots_read_pair_numbers = out_prefix + "_read_pair_number_boxplot.pdf"

# PDF file with boxplots for distributions of interaction distances
pdf_name_boxplots_interaction_distances = out_prefix + "_interaction_distance_boxplot.pdf"


### 1st pass: Collect information about DI, UIE and UII
#######################################################

print("[INFO] 1st pass: Collect information about DI, UIE and UII ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Collect read pair numbers for DI, UIE and UII
        if interaction_category == 'DI':
            dir_inter_num += 1
            if enrichment_pair_tag == 'AA':
                dir_inter_aa_num +=1
                dir_inter_aa_rp_array.append(rp_total)
                dir_inter_aa_dist_array.append(i_dist)
                if rp_total not in rp_dict_aa:
                    rp_dict_aa[rp_total] = 1
                else:
                    rp_dict_aa[rp_total] +=1
                # Write directed interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, 'DIAA') + "\n")
            elif enrichment_pair_tag == 'AI' or  enrichment_pair_tag == 'IA':
                dir_inter_ai_num += 1
                dir_inter_ai_rp_array.append(rp_total)
                dir_inter_ai_dist_array.append(i_dist)
                if rp_total not in rp_dict_ai:
                    rp_dict_ai[rp_total] = 1
                else:
                    rp_dict_ai[rp_total] +=1
                # Write directed interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, 'DIAI') + "\n")
            elif enrichment_pair_tag == 'II':
                dir_inter_ii_num += 1
                dir_inter_ii_rp_array.append(rp_total)
                dir_inter_ii_dist_array.append(i_dist)
                if rp_total not in rp_dict_ii:
                    rp_dict_ii[rp_total] = 1
                else:
                    rp_dict_ii[rp_total] +=1
                # Write directed interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, 'DIII') + "\n")

        elif interaction_category == 'UIE' or interaction_category == 'UII':
            undir_inter_num += 1
            if enrichment_pair_tag == 'AA':
                undir_inter_aa_num +=1
                undir_inter_aa_rp_array.append(rp_total)
                undir_inter_aa_dist_array.append(i_dist)
            elif enrichment_pair_tag == 'AI' or  enrichment_pair_tag == 'IA':
                undir_inter_ai_num += 1
                undir_inter_ai_rp_array.append(rp_total)
                undir_inter_ai_dist_array.append(i_dist)
            elif enrichment_pair_tag == 'II':
                undir_inter_ii_num += 1
                undir_inter_ii_rp_array.append(rp_total)
                undir_inter_ii_dist_array.append(i_dist)
        else:
            print("[Error] Interaction category must be either \'DI\', \'UIE\' or \'UII\'!")
            exit(1)

        line = fp.readline()


### Get Q1 and Q3 for AA, AI and II
###################################

# Read pair numbers
dir_inter_aa_rp_q1 = int(numpy.quantile(dir_inter_aa_rp_array, 0.25))
dir_inter_ai_rp_q1 = int(numpy.quantile(dir_inter_ai_rp_array, 0.25))
dir_inter_ii_rp_q1 = int(numpy.quantile(dir_inter_ii_rp_array, 0.25))

dir_inter_aa_rp_q3 = int(numpy.quantile(dir_inter_aa_rp_array, 0.75))
dir_inter_ai_rp_q3 = int(numpy.quantile(dir_inter_ai_rp_array, 0.75))
dir_inter_ii_rp_q3 = int(numpy.quantile(dir_inter_ii_rp_array, 0.75))

# Interaction distances
dir_inter_aa_dist_q1 = int(numpy.quantile(dir_inter_aa_dist_array, 0.25))
dir_inter_ai_dist_q1 = int(numpy.quantile(dir_inter_ai_dist_array, 0.25))
dir_inter_ii_dist_q1 = int(numpy.quantile(dir_inter_ii_dist_array, 0.25))

dir_inter_aa_dist_q3 = int(numpy.quantile(dir_inter_aa_dist_array, 0.75))
dir_inter_ai_dist_q3 = int(numpy.quantile(dir_inter_ai_dist_array, 0.75))
dir_inter_ii_dist_q3 = int(numpy.quantile(dir_inter_ii_dist_array, 0.75))


### 2nd pass: Select undirected reference interactions
######################################################

print("[INFO] 2nd pass: Select undirected reference interactions ...")

print("\t[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Report progress
        n_interaction_total += 1
        if n_interaction_total % 1000000 == 0:
            print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Adjust reference interactions to interaction size
        if adjust_to_interaction_distance:
            i_dist_range_aa_ok = dir_inter_aa_dist_q1 <= i_dist and i_dist <= dir_inter_aa_dist_q3
            i_dist_range_ai_ok = dir_inter_ai_dist_q1 <= i_dist and i_dist <= dir_inter_ai_dist_q3
            i_dist_range_ii_ok = dir_inter_ii_dist_q1 <= i_dist and i_dist <= dir_inter_ii_dist_q3
        else:
            i_dist_range_aa_ok = True
            i_dist_range_ai_ok = True
            i_dist_range_ii_ok = True

        # Select reference interactions (1) using Q1 and Q2 for read pair numbers in DI
        if interaction_category == 'UIE' or interaction_category == 'UII':
            interaction_category_tag = 'NA'
            if enrichment_pair_tag == 'AA' and dir_inter_aa_rp_q1 < rp_total and rp_total < dir_inter_aa_rp_q3 and i_dist_range_aa_ok:
                undir_ref_1_inter_aa_rp_array.append(rp_total)
                undir_ref_1_inter_aa_dist_array.append(i_dist)
                undir_ref_1_inter_num +=1
                undir_ref_1_inter_aa_num += 1
                interaction_category_tag = 'UIRAA'
            elif (enrichment_pair_tag == 'AI' or enrichment_pair_tag == 'IA') and dir_inter_ai_rp_q1 < rp_total and rp_total < dir_inter_ai_rp_q3 and i_dist_range_ai_ok:
                undir_ref_1_inter_ai_rp_array.append(rp_total)
                undir_ref_1_inter_ai_dist_array.append(i_dist)
                undir_ref_1_inter_num += 1
                undir_ref_1_inter_ai_num += 1
                interaction_category_tag = 'UIRAI'
            elif enrichment_pair_tag == 'II' and dir_inter_ii_rp_q1 < rp_total and rp_total < dir_inter_ii_rp_q3 and i_dist_range_ii_ok:
                undir_ref_1_inter_ii_rp_array.append(rp_total)
                undir_ref_1_inter_ii_dist_array.append(i_dist)
                undir_ref_1_inter_num += 1
                undir_ref_1_inter_ii_num += 1
                interaction_category_tag = 'UIRII'

            if selection_approach == 'q13' and interaction_category_tag != 'NA':
                # Override interaction category tag in column 3 and write interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, interaction_category_tag) + "\n")

        # Select reference interactions (2) using dictionaries for read pair numbers DI
        if interaction_category == 'UIE' or interaction_category == 'UII':
            interaction_category_tag = 'NA'
            if enrichment_pair_tag == 'AA' and rp_total in rp_dict_aa and 0 < rp_dict_aa[rp_total] and i_dist_range_aa_ok:
                rp_dict_aa[rp_total] -= 1
                undir_ref_2_inter_aa_rp_array.append(rp_total)
                undir_ref_2_inter_aa_dist_array.append(i_dist)
                undir_ref_2_inter_num +=1
                undir_ref_2_inter_aa_num += 1
                interaction_category_tag = 'UIRAA'
            elif (enrichment_pair_tag == 'AI' or enrichment_pair_tag == 'IA') and rp_total in rp_dict_ai and 0 < rp_dict_ai[rp_total] and i_dist_range_ai_ok:
                rp_dict_ai[rp_total] -= 1
                undir_ref_2_inter_ai_rp_array.append(rp_total)
                undir_ref_2_inter_ai_dist_array.append(i_dist)
                undir_ref_2_inter_num += 1
                undir_ref_2_inter_ai_num += 1
                interaction_category_tag = 'UIRAI'
            elif enrichment_pair_tag == 'II' and rp_total in rp_dict_ii and 0 < rp_dict_ii[rp_total] and i_dist_range_ii_ok:
                rp_dict_ii[rp_total] -= 1
                undir_ref_2_inter_ii_rp_array.append(rp_total)
                undir_ref_2_inter_ii_dist_array.append(i_dist)
                undir_ref_2_inter_num += 1
                undir_ref_2_inter_ii_num += 1
                interaction_category_tag = 'UIRII'

            if selection_approach == 'exact' and interaction_category_tag != 'NA':
                # Override interaction category tag in column 3 and write interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, interaction_category_tag) + "\n")

        line = fp.readline()

enhanced_interaction_stream_output.close()


### Output some statistics
##########################

dir_inter_aa_rp_median = int(numpy.quantile(dir_inter_aa_rp_array, 0.50))
dir_inter_ai_rp_median = int(numpy.quantile(dir_inter_ai_rp_array, 0.50))
dir_inter_ii_rp_median = int(numpy.quantile(dir_inter_ii_rp_array, 0.50))

dir_inter_aa_percentage = "{0:.2f}".format(100 * dir_inter_aa_num / dir_inter_num)
dir_inter_ai_percentage = "{0:.2f}".format(100 * dir_inter_ai_num / dir_inter_num)
dir_inter_ii_percentage = "{0:.2f}".format(100 * dir_inter_ii_num / dir_inter_num)

undir_inter_aa_rp_median = int(numpy.quantile(undir_inter_aa_rp_array, 0.50))
undir_inter_ai_rp_median = int(numpy.quantile(undir_inter_ai_rp_array, 0.50))
undir_inter_ii_rp_median = int(numpy.quantile(undir_inter_ii_rp_array, 0.50))

undir_inter_aa_percentage = "{0:.2f}".format(100 * undir_inter_aa_num / undir_inter_num)
undir_inter_ai_percentage = "{0:.2f}".format(100 * undir_inter_ai_num / undir_inter_num)
undir_inter_ii_percentage = "{0:.2f}".format(100 * undir_inter_ii_num / undir_inter_num)

undir_ref_1_inter_aa_rp_median = int(numpy.quantile(undir_ref_1_inter_aa_rp_array, 0.50))
undir_ref_1_inter_ai_rp_median = int(numpy.quantile(undir_ref_1_inter_ai_rp_array, 0.50))
undir_ref_1_inter_ii_rp_median = int(numpy.quantile(undir_ref_1_inter_ii_rp_array, 0.50))

undir_ref_1_inter_aa_percentage = "{0:.2f}".format(100 * undir_ref_1_inter_aa_num / undir_ref_1_inter_num)
undir_ref_1_inter_ai_percentage = "{0:.2f}".format(100 * undir_ref_1_inter_ai_num / undir_ref_1_inter_num)
undir_ref_1_inter_ii_percentage = "{0:.2f}".format(100 * undir_ref_1_inter_ii_num / undir_ref_1_inter_num)

undir_ref_2_inter_aa_rp_median = int(numpy.quantile(undir_ref_2_inter_aa_rp_array, 0.50))
undir_ref_2_inter_ai_rp_median = int(numpy.quantile(undir_ref_2_inter_ai_rp_array, 0.50))
undir_ref_2_inter_ii_rp_median = int(numpy.quantile(undir_ref_2_inter_ii_rp_array, 0.50))

undir_ref_2_inter_aa_percentage = "{0:.2f}".format(100 * undir_ref_2_inter_aa_num / undir_ref_2_inter_num)
undir_ref_2_inter_ai_percentage = "{0:.2f}".format(100 * undir_ref_2_inter_ai_num / undir_ref_2_inter_num)
undir_ref_2_inter_ii_percentage = "{0:.2f}".format(100 * undir_ref_2_inter_ii_num / undir_ref_2_inter_num)

print()
print("Total number of directed interactions: " + str(dir_inter_num))
print("\tWithin 'AA': " + str(dir_inter_aa_num) + " (" + str(dir_inter_aa_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(dir_inter_aa_rp_median))
print("\tWithin 'AI': " + str(dir_inter_ai_num) + " (" + str(dir_inter_ai_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(dir_inter_ai_rp_median))
print("\tWithin 'II': " + str(dir_inter_ii_num) + " (" + str(dir_inter_ii_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(dir_inter_ii_rp_median))

print()
print("Number of exclusive and inclusive undirected interactions: " + str(undir_inter_num))
print("\tWithin 'AA': " + str(undir_inter_aa_num) + " (" + str(undir_inter_aa_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_inter_aa_rp_median))
print("\tWithin 'AI': " + str(undir_inter_ai_num) + " (" + str(undir_inter_ai_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_inter_ai_rp_median))
print("\tWithin 'II': " + str(undir_inter_ii_num) + " (" + str(undir_inter_ii_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_inter_ii_rp_median))

print()
print("Number of exclusive and inclusive undirected reference (1) interactions: " + str(undir_ref_1_inter_num))
print("\tWithin 'AA': " + str(undir_ref_1_inter_aa_num) + " (" + str(undir_ref_1_inter_aa_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_1_inter_aa_rp_median))
print("\tWithin 'AI': " + str(undir_ref_1_inter_ai_num) + " (" + str(undir_ref_1_inter_ai_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_1_inter_ai_rp_median))
print("\tWithin 'II': " + str(undir_ref_1_inter_ii_num) + " (" + str(undir_ref_1_inter_ii_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_1_inter_ii_rp_median))

print()
print("Number of exclusive and inclusive undirected reference (2) interactions: " + str(undir_ref_2_inter_num))
print("\tWithin 'AA': " + str(undir_ref_2_inter_aa_num) + " (" + str(undir_ref_2_inter_aa_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_2_inter_aa_rp_median))
print("\tWithin 'AI': " + str(undir_ref_2_inter_ai_num) + " (" + str(undir_ref_2_inter_ai_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_2_inter_ai_rp_median))
print("\tWithin 'II': " + str(undir_ref_2_inter_ii_num) + " (" + str(undir_ref_2_inter_ii_percentage) + "%)")
print("\t\tMedian number of read pairs: " + str(undir_ref_2_inter_ii_rp_median))

print()


### Create boxplots for read pair numbers
#########################################

data = [
    dir_inter_aa_rp_array, dir_inter_ai_rp_array, dir_inter_ii_rp_array,
    undir_inter_aa_rp_array,undir_inter_ai_rp_array, undir_inter_ii_rp_array,
    undir_ref_1_inter_aa_rp_array, undir_ref_1_inter_ai_rp_array, undir_ref_1_inter_ii_rp_array,
    undir_ref_2_inter_aa_rp_array, undir_ref_2_inter_ai_rp_array, undir_ref_2_inter_ii_rp_array
]
labels = ['AA', 'AI', 'II', 'AA', 'AI', 'II','AA', 'AI', 'II','AA', 'AI', 'II']
fig1, ax1 = plt.subplots()
ax1.set_title('Distributions of read pair numbers for DI, U, UR1 and UR2 (' + out_prefix + ')')
ax1.boxplot(data, showfliers=False, labels=labels)

box = ax1.boxplot(data, showfliers=False, labels=labels, patch_artist=True)

colors = ['blue', 'blue', 'blue', 'gray', 'gray', 'gray', 'lavender', 'lavender', 'lavender', 'pink', 'pink', 'pink']

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

lab = "DI - AA: " + str(dir_inter_aa_num) + " (" + dir_inter_aa_percentage + "%), AI: " + str(dir_inter_ai_num) + " (" + dir_inter_ai_percentage + "%), II: " + str(dir_inter_ii_num) + " (" + dir_inter_ii_percentage + "%)"
di_patch = mpatches.Patch(color='blue', label=lab)
lab = "U - AA: " + str(undir_inter_aa_num) + " (" + undir_inter_aa_percentage + "%), AI: " + str(undir_inter_ai_num) + " (" + undir_inter_ai_percentage + "%), II: " + str(undir_inter_ii_num) + " (" + undir_inter_ii_percentage + "%)"
u_patch = mpatches.Patch(color='gray', label=lab)
lab = "UR 1 - AA: " + str(undir_ref_1_inter_aa_num) + " (" + undir_ref_1_inter_aa_percentage + "%), AI: " + str(undir_ref_1_inter_ai_num) + " (" + undir_ref_1_inter_ai_percentage + "%), II: " + str(undir_ref_1_inter_ii_num) + " (" + undir_ref_1_inter_ii_percentage + "%)"
ur_patch = mpatches.Patch(color='lavender', label=lab)
lab = "UR 2 - AA: " + str(undir_ref_2_inter_aa_num) + " (" + undir_ref_2_inter_aa_percentage + "%), AI: " + str(undir_ref_2_inter_ai_num) + " (" + undir_ref_2_inter_ai_percentage + "%), II: " + str(undir_ref_2_inter_ii_num) + " (" + undir_ref_2_inter_ii_percentage + "%)"
ur2_patch = mpatches.Patch(color='pink', label=lab)
plt.legend(handles=[di_patch, u_patch, ur_patch, ur2_patch])
ax1.set_xlabel('Enrichment state pair tag')
ax1.set_ylabel('Read pair number')
plt.grid(True)
fig1.set_size_inches(10,5)
plt.savefig(pdf_name_boxplots_read_pair_numbers)
plt.close()


### Create boxplots for interaction distances
#############################################

data = [
    dir_inter_aa_dist_array, dir_inter_ai_dist_array, dir_inter_ii_dist_array,
    undir_inter_aa_dist_array,undir_inter_ai_dist_array, undir_inter_ii_dist_array,
    undir_ref_1_inter_aa_dist_array, undir_ref_1_inter_ai_dist_array, undir_ref_1_inter_ii_dist_array,
    undir_ref_2_inter_aa_dist_array, undir_ref_2_inter_ai_dist_array, undir_ref_2_inter_ii_dist_array
]
labels = ['AA', 'AI', 'II', 'AA', 'AI', 'II','AA', 'AI', 'II','AA', 'AI', 'II']
fig1, ax1 = plt.subplots()
ax1.set_title('Distributions of interaction distances for DI, U, UR 1 and UR 2 (' + out_prefix + ')')
ax1.boxplot(data, showfliers=False, labels=labels)

box = ax1.boxplot(data, showfliers=False, labels=labels, patch_artist=True)

colors = ['blue', 'blue', 'blue', 'gray', 'gray', 'gray', 'lavender', 'lavender', 'lavender', 'pink', 'pink', 'pink']

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

lab = "DI - AA: " + str(dir_inter_aa_num) + " (" + dir_inter_aa_percentage + "%), AI: " + str(dir_inter_ai_num) + " (" + dir_inter_ai_percentage + "%), II: " + str(dir_inter_ii_num) + " (" + dir_inter_ii_percentage + "%)"
di_patch = mpatches.Patch(color='blue', label=lab)
lab = "U - AA: " + str(undir_inter_aa_num) + " (" + undir_inter_aa_percentage + "%), AI: " + str(undir_inter_ai_num) + " (" + undir_inter_ai_percentage + "%), II: " + str(undir_inter_ii_num) + " (" + undir_inter_ii_percentage + "%)"
u_patch = mpatches.Patch(color='gray', label=lab)
lab = "UR 1 - AA: " + str(undir_ref_1_inter_aa_num) + " (" + undir_ref_1_inter_aa_percentage + "%), AI: " + str(undir_ref_1_inter_ai_num) + " (" + undir_ref_1_inter_ai_percentage + "%), II: " + str(undir_ref_1_inter_ii_num) + " (" + undir_ref_1_inter_ii_percentage + "%)"
ur_patch = mpatches.Patch(color='lavender', label=lab)
lab = "UR 2 - AA: " + str(undir_ref_2_inter_aa_num) + " (" + undir_ref_2_inter_aa_percentage + "%), AI: " + str(undir_ref_2_inter_ai_num) + " (" + undir_ref_2_inter_ai_percentage + "%), II: " + str(undir_ref_2_inter_ii_num) + " (" + undir_ref_2_inter_ii_percentage + "%)"
ur2_patch = mpatches.Patch(color='pink', label=lab)
plt.legend(handles=[di_patch, u_patch, ur_patch, ur2_patch])
ax1.set_xlabel('Enrichment state pair tag')
ax1.set_ylabel('Interaction distance')
plt.grid(True)
fig1.set_size_inches(10,5)
plt.savefig(pdf_name_boxplots_interaction_distances)
plt.close()
