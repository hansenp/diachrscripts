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
and another one for interaction distances:

   '<OUT_PREFIX>_read_pair_number_boxplot.pdf'
   '<OUT_PREFIX>_interaction_distance_boxplot.pdf'

Furthermore, barplots for the proportions of interactions and associated digests are created:

   '<OUT_PREFIX>_interaction_enrichment_pair_tags_barplot.pdf'
   '<OUT_PREFIX>__digest_enrichment_tags_barplot.pdf'

"""


import argparse
import gzip
import diachrscripts_toolkit
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.backends.backend_pdf


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

# Second quartiles of read pair number distribution within AA, AI and II
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

# Second quartiles of read pair number distribution within AA, AI and II
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

# Sets for digests involved directed and undirected interactions
dir_a_dig_set = set()
dir_i_dig_set = set()
undir_a_dig_set = set()
undir_i_dig_set = set()
undir_ref_1_a_dig_set = set()
undir_ref_1_i_dig_set = set()
undir_ref_2_a_dig_set = set()
undir_ref_2_i_dig_set = set()

# Output files
##############

# Enhanced interaction file with interaction categories: DI, UIRAA, UIRAI and UIRII
enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_with_di_and_uir.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')

# PDF file with boxplots for distributions of read pair numbers
pdf_name_boxplots_read_pair_numbers = out_prefix + "_read_pair_number_boxplot.pdf"

# PDF file with boxplots for distributions of interaction distances
pdf_name_boxplots_interaction_distances = out_prefix + "_interaction_distance_boxplot.pdf"

# PDF file with three barplots showing the proportions of interactions within the enrichment pair categories AA, AI and II
pdf_name_barplots_interaction_enrichment_pair_tags = out_prefix + "_interaction_enrichment_pair_tags_barplot.pdf"

# Prepare stream for warnings
tab_file_warnings_output = out_prefix + "_warnings_di_and_uri.tsv"
tab_stream_warnings_output = open(tab_file_warnings_output, 'wt')

# Prepare stream for output of statistics
tab_file_stats_output = out_prefix + "_stats_di_and_uri.tsv"
tab_stream_stats_output = open(tab_file_stats_output, 'wt')



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

        enrichment_tag_dig_1 = enrichment_pair_tag[0]
        enrichment_tag_dig_2 = enrichment_pair_tag[1]

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

            if enrichment_tag_dig_1 == 'A':
                dir_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            else:
                dir_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            if enrichment_tag_dig_2 == 'A':
                dir_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            else:
                dir_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

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
            if enrichment_tag_dig_1 == 'A':
                undir_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            else:
                undir_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
            if enrichment_tag_dig_2 == 'A':
                undir_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            else:
                undir_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
        else:
            print("[Error] Interaction category must be either \'DI\', \'UIE\' or \'UII\'!")
            exit(1)

        line = fp.readline()

fp.close()
print("\t[INFO] done ...")

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

        enrichment_tag_dig_1 = enrichment_pair_tag[0]
        enrichment_tag_dig_2 = enrichment_pair_tag[1]

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
                undir_ref_1_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                undir_ref_1_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            elif (enrichment_pair_tag == 'AI' or enrichment_pair_tag == 'IA') and dir_inter_ai_rp_q1 < rp_total and rp_total < dir_inter_ai_rp_q3 and i_dist_range_ai_ok:
                undir_ref_1_inter_ai_rp_array.append(rp_total)
                undir_ref_1_inter_ai_dist_array.append(i_dist)
                undir_ref_1_inter_num += 1
                undir_ref_1_inter_ai_num += 1
                interaction_category_tag = 'UIRAI'
                if enrichment_tag_dig_1 == 'A':
                    undir_ref_1_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                else:
                    undir_ref_1_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                if enrichment_tag_dig_2 == 'A':
                    undir_ref_1_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
                else:
                    undir_ref_1_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            elif enrichment_pair_tag == 'II' and dir_inter_ii_rp_q1 < rp_total and rp_total < dir_inter_ii_rp_q3 and i_dist_range_ii_ok:
                undir_ref_1_inter_ii_rp_array.append(rp_total)
                undir_ref_1_inter_ii_dist_array.append(i_dist)
                undir_ref_1_inter_num += 1
                undir_ref_1_inter_ii_num += 1
                interaction_category_tag = 'UIRII'
                undir_ref_1_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                undir_ref_1_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

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
                undir_ref_2_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                undir_ref_2_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            elif (enrichment_pair_tag == 'AI' or enrichment_pair_tag == 'IA') and rp_total in rp_dict_ai and 0 < rp_dict_ai[rp_total] and i_dist_range_ai_ok:
                rp_dict_ai[rp_total] -= 1
                undir_ref_2_inter_ai_rp_array.append(rp_total)
                undir_ref_2_inter_ai_dist_array.append(i_dist)
                undir_ref_2_inter_num += 1
                undir_ref_2_inter_ai_num += 1
                interaction_category_tag = 'UIRAI'
                if enrichment_tag_dig_1 == 'A':
                    undir_ref_2_a_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                else:
                    undir_ref_2_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                if enrichment_tag_dig_2 == 'A':
                    undir_ref_2_a_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
                else:
                    undir_ref_2_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
            elif enrichment_pair_tag == 'II' and rp_total in rp_dict_ii and 0 < rp_dict_ii[rp_total] and i_dist_range_ii_ok:
                rp_dict_ii[rp_total] -= 1
                undir_ref_2_inter_ii_rp_array.append(rp_total)
                undir_ref_2_inter_ii_dist_array.append(i_dist)
                undir_ref_2_inter_num += 1
                undir_ref_2_inter_ii_num += 1
                interaction_category_tag = 'UIRII'
                undir_ref_2_i_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                undir_ref_2_i_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            if selection_approach == 'exact' and interaction_category_tag != 'NA':
                # Override interaction category tag in column 3 and write interaction to file
                enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, interaction_category_tag) + "\n")

        line = fp.readline()

enhanced_interaction_stream_output.close()
fp.close()
print("\t[INFO] done ...")

print("[INFO] Checking whether reference interactions are missing for some n")

#print("\t[INFO] AA interactions")
tab_stream_warnings_output.write("AA_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_aa = 0
for x in rp_dict_aa:
    if 0 < rp_dict_aa[x]:
        n_missing_aa += rp_dict_aa[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(rp_dict_aa[x]) + "\n")
        #print("\t\t[WARNING] For " + str(rp_dict_aa[x]) + " directed interactions with " + str(x) + " read pairs no undirected reference interaction could be selected.")

#print("\t[INFO] AI interactions")
tab_stream_warnings_output.write("AI_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_ai = 0
for x in rp_dict_ai:
    if 0 < rp_dict_ai[x]:
        n_missing_ai += rp_dict_ai[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(rp_dict_ai[x]) + "\n")
        #print("\t\t[WARNING] For " + str(rp_dict_ai[x]) + " directed interactions with " + str(x) + " read pairs no undirected reference interaction could be selected.")

#print("\t[INFO] II interactions")
tab_stream_warnings_output.write("II_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_ii = 0
for x in rp_dict_ii:
    if 0 < rp_dict_ii[x]:
        n_missing_ii += rp_dict_ii[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(rp_dict_ii[x]) + "\n")
        #print("\t\t[WARNING] For " + str(rp_dict_ii[x]) + " directed interactions with " + str(x) + " read pairs no undirected reference interaction could be selected.")

print("\t[INFO] Summary for AA, AI and II")
print("\t\t[WARNING] For " + str(n_missing_aa) + " out of " + str(dir_inter_aa_num) + " directed AA interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_ai) + " out of " + str(dir_inter_ai_num) + " directed AI interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_ii) + " out of " + str(dir_inter_ii_num) + " directed II interactions no undirected reference interaction could be selected.")

tab_stream_warnings_output.close()
print("\t[INFO] done ...")

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

dir_a_dig_num = len(dir_a_dig_set)                   # Number of active digests involved in directed interactions
dir_i_dig_num = len(dir_i_dig_set)                   # Number of inactive digests involved in directed interactions
percentage_dir_a_dig = 100 * dir_a_dig_num / (dir_a_dig_num + dir_i_dig_num)

undir_a_dig_num = len(undir_a_dig_set)               # Number of active digests involved in undirected interactions
undir_i_dig_num = len(undir_i_dig_set)               # Number of inactive digests involved in undirected interactions
percentage_undir_a_dig = 100 * undir_a_dig_num / (undir_a_dig_num + undir_i_dig_num)

undir_ref_1_a_dig_num = len(undir_ref_1_a_dig_set)   # Number of active digests involved in undirected interactions
undir_ref_1_i_dig_num = len(undir_ref_1_i_dig_set)   # Number of inactive digests involved in undirected interactions
percentage_undir_ref_1_a_dig = 100 * undir_ref_1_a_dig_num / (undir_ref_1_a_dig_num + undir_ref_1_i_dig_num)

undir_ref_2_a_dig_num = len(undir_ref_2_a_dig_set)   # Number of active digests involved in undirected interactions
undir_ref_2_i_dig_num = len(undir_ref_2_i_dig_set)   # Number of inactive digests involved in undirected interactions
percentage_undir_ref_2_a_dig = 100 * undir_ref_2_a_dig_num / (undir_ref_2_a_dig_num + undir_ref_2_i_dig_num)

print()
print("Number of active (A) digests involved in DI: " + str(dir_a_dig_num) + " (" + "{0:.2f}".format(percentage_dir_a_dig) + "%)")
print("Number of inactive (I) digests involved in DI: " + str(dir_i_dig_num))
print()
print("Number of active (A) digests involved in U: " + str(undir_a_dig_num) + " (" + "{0:.2f}".format(percentage_undir_a_dig) + "%)")
print("Number of inactive (I) digests involved in U: " + str(undir_i_dig_num))
print()
print("Number of active (A) digests involved in UR1: " + str(undir_ref_1_a_dig_num) + " (" + "{0:.2f}".format(percentage_undir_ref_1_a_dig) + "%)")
print("Number of inactive (I) digests involved in UR1: " + str(undir_ref_1_i_dig_num))
print()
print("Number of active (A) digests involved in UR2: " + str(undir_ref_2_a_dig_num) + " (" + "{0:.2f}".format(percentage_undir_ref_2_a_dig) + "%)")
print("Number of inactive (I) digests involved in UR2: " + str(undir_ref_2_i_dig_num))
print()

dir_inter_aa_dist_median = int(numpy.quantile(dir_inter_aa_dist_array, 0.50))
dir_inter_ai_dist_median = int(numpy.quantile(dir_inter_ai_dist_array, 0.50))
dir_inter_ii_dist_median = int(numpy.quantile(dir_inter_ii_dist_array, 0.50))

print("Median distances of directed interactions:")
print("\tAA: " + str(dir_inter_aa_dist_median))
print("\tAI: " + str(dir_inter_ai_dist_median))
print("\tII: " + str(dir_inter_ii_dist_median))

undir_inter_aa_dist_median = int(numpy.quantile(undir_inter_aa_dist_array, 0.50))
undir_inter_ai_dist_median = int(numpy.quantile(undir_inter_ai_dist_array, 0.50))
undir_inter_ii_dist_median = int(numpy.quantile(undir_inter_ii_dist_array, 0.50))

print("Median distances of undirected interactions:")
print("\tAA: " + str(undir_inter_aa_dist_median))
print("\tAI: " + str(undir_inter_ai_dist_median))
print("\tII: " + str(undir_inter_ii_dist_median))

undir_ref_1_inter_aa_dist_median = int(numpy.quantile(undir_ref_1_inter_aa_dist_array, 0.50))
undir_ref_1_inter_ai_dist_median = int(numpy.quantile(undir_ref_1_inter_ai_dist_array, 0.50))
undir_ref_1_inter_ii_dist_median = int(numpy.quantile(undir_ref_1_inter_ii_dist_array, 0.50))

print("Median distances of undirected reference interactions (q13):")
print("\tAA: " + str(undir_ref_1_inter_aa_dist_median))
print("\tAI: " + str(undir_ref_1_inter_ai_dist_median))
print("\tII: " + str(undir_ref_1_inter_ii_dist_median))

undir_ref_2_inter_aa_dist_median = int(numpy.quantile(undir_ref_2_inter_aa_dist_array, 0.50))
undir_ref_2_inter_ai_dist_median = int(numpy.quantile(undir_ref_2_inter_ai_dist_array, 0.50))
undir_ref_2_inter_ii_dist_median = int(numpy.quantile(undir_ref_2_inter_ii_dist_array, 0.50))

print("Median distances of undirected reference interactions (exact):")
print("\tAA: " + str(undir_ref_2_inter_aa_dist_median))
print("\tAI: " + str(undir_ref_2_inter_ai_dist_median))
print("\tII: " + str(undir_ref_2_inter_ii_dist_median))

print()

tab_stream_stats_output.write(

    "out_prefix" + "\t" +                               # Prefix for output

    "dir_inter_num" + "\t" +                            # Total number of directed interactions

    "dir_inter_aa_num" + "\t" +                         # Number of directed interactions within AA
    "dir_inter_ai_num" + "\t" +                         # Number of directed interactions within AI
    "dir_inter_ii_num" + "\t" +                         # Number of directed interactions within II

    "dir_inter_aa_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within AA
    "dir_inter_ai_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within AI
    "dir_inter_ii_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within II

    "undir_inter_num" + "\t" +                          # Total number of undirected interactions

    "undir_inter_aa_num" + "\t" +                       # Number of undirected interactions within AA
    "undir_inter_ai_num" + "\t" +                       # Number of undirected interactions within AI
    "undir_inter_ii_num" + "\t" +                       # Number of undirected interactions within II

    "undir_inter_aa_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within AA
    "undir_inter_ai_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within AI
    "undir_inter_ii_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within II

    "undir_ref_1_inter_num" + "\t" +                    # Total number of undirected reference interactions (q13)

    "undir_ref_1_inter_aa_num" + "\t" +                 # Number of undirected reference interactions within AA
    "undir_ref_1_inter_ai_num" + "\t" +                 # Number of undirected reference interactions within AI
    "undir_ref_1_inter_ii_num" + "\t" +                 # Number of undirected reference interactions within II

    "undir_ref_1_inter_aa_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AA
    "undir_ref_1_inter_ai_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AI
    "undir_ref_1_inter_ii_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within II

    "undir_ref_2_inter_num" + "\t" +                    # Total number of undirected interactions (exact)

    "undir_ref_2_inter_aa_num" + "\t" +                 # Number of undirected reference interactions within AA
    "undir_ref_2_inter_ai_num" + "\t" +                 # Number of undirected reference interactions within AI
    "undir_ref_2_inter_ii_num" + "\t" +                 # Number of undirected reference interactions within II

    "n_missing_aa" + "\t" +                             # Number of missing reference interactions within AA
    "n_missing_ai" + "\t" +                             # Number of missing reference interactions within AI
    "n_missing_ii" + "\t" +                             # Number of missing reference interactions within II

    "undir_ref_2_inter_aa_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AA
    "undir_ref_2_inter_ai_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AI
    "undir_ref_2_inter_ii_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within II

    "dir_a_dig_num" + "\t" +                            # Number of active digests involved in directed interactions
    "dir_i_dig_num" + "\t" +                            # Number of inactive digests involved in directed interactions

    "undir_a_dig_num" + "\t" +                          # Number of active digests involved in undirected interactions
    "undir_i_dig_num" + "\t" +                          # Number of inactive digests involved in undirected interactions

    "undir_ref_1_a_dig_num" + "\t" +                    # Number of active digests involved in undirected reference interactions (q13)
    "undir_ref_1_i_dig_num" + "\t" +                    # Number of inactive digests involved in undirected reference interactions (q13)

    "undir_ref_2_a_dig_num" + "\t" +                    # Number of active digests involved in undirected reference interactions (exact)
    "undir_ref_2_i_dig_num" + "\t" +                    # Number of inactive digests involved in undirected reference interactions (exact)

    "dir_inter_aa_dist_median" + "\t" +                 # Median interaction distance of directed interactions within AA
    "dir_inter_ai_dist_median" + "\t" +                 # Median interaction distance of directed interactions within AI
    "dir_inter_ii_dist_median" + "\t" +                 # Median interaction distance of directed interactions within II

    "undir_inter_aa_dist_median" + "\t" +               # Median interaction distance of undirected interactions within AA
    "undir_inter_ai_dist_median" + "\t" +               # Median interaction distance of undirected interactions within AI
    "undir_inter_ii_dist_median" + "\t" +               # Median interaction distance of undirected interactions within II

    "undir_ref_1_inter_aa_dist_median" + "\t" +         # Median interaction distance of undirected reference interactions within AA (q13)
    "undir_ref_1_inter_ai_dist_median" + "\t" +         # Median interaction distance of undirected reference interactions within AI (q13)
    "undir_ref_1_inter_ii_dist_median" + "\t" +         # Median interaction distance of undirected reference interactions within II (q13)

    "undir_ref_2_inter_aa_dist_median" + "\t" +         # Median interaction distance of undirected reference interactions within AA (exact)
    "undir_ref_2_inter_ai_dist_median" + "\t" +         # Median interaction distance of undirected reference interactions within AI (exact)
    "undir_ref_2_inter_ii_dist_median" +                # Median interaction distance of undirected reference interactions within II (exact)

    "\n"
)
tab_stream_stats_output.write(

    str(out_prefix) + "\t" +

    str(dir_inter_num) + "\t" +

    str(dir_inter_aa_num) + "\t" +
    str(dir_inter_ai_num) + "\t" +
    str(dir_inter_ii_num) + "\t" +

    str(dir_inter_aa_rp_median) + "\t" +
    str(dir_inter_ai_rp_median) + "\t" +
    str(dir_inter_ii_rp_median) + "\t" +

    str(undir_inter_num) + "\t" +

    str(undir_inter_aa_num) + "\t" +
    str(undir_inter_ai_num) + "\t" +
    str(undir_inter_ii_num) + "\t" +

    str(undir_inter_aa_rp_median) + "\t" +
    str(undir_inter_ai_rp_median) + "\t" +
    str(undir_inter_ii_rp_median) + "\t" +

    str(undir_ref_1_inter_num) + "\t" +

    str(undir_ref_1_inter_aa_num) + "\t" +
    str(undir_ref_1_inter_ai_num) + "\t" +
    str(undir_ref_1_inter_ii_num) + "\t" +

    str(undir_ref_1_inter_aa_rp_median) + "\t" +
    str(undir_ref_1_inter_ai_rp_median) + "\t" +
    str(undir_ref_1_inter_ii_rp_median) + "\t" +

    str(undir_ref_2_inter_num) + "\t" +

    str(undir_ref_2_inter_aa_num) + "\t" +
    str(undir_ref_2_inter_ai_num) + "\t" +
    str(undir_ref_2_inter_ii_num) + "\t" +

    str(n_missing_aa) + "\t" +
    str(n_missing_ai) + "\t" +
    str(n_missing_ii) + "\t" +

    str(undir_ref_2_inter_aa_rp_median) + "\t" +
    str(undir_ref_2_inter_ai_rp_median) + "\t" +
    str(undir_ref_2_inter_ii_rp_median) + "\t" +

    str(dir_a_dig_num) + "\t" +
    str(dir_i_dig_num) + "\t" +

    str(undir_a_dig_num) + "\t" +
    str(undir_i_dig_num) + "\t" +

    str(undir_ref_1_a_dig_num) + "\t" +
    str(undir_ref_1_i_dig_num) + "\t" +

    str(undir_ref_2_a_dig_num) + "\t" +
    str(undir_ref_2_i_dig_num) + "\t" +

    str(dir_inter_aa_dist_median) + "\t" +
    str(dir_inter_ai_dist_median) + "\t" +
    str(dir_inter_ii_dist_median) + "\t" +

    str(undir_inter_aa_dist_median) + "\t" +
    str(undir_inter_ai_dist_median) + "\t" +
    str(undir_inter_ii_dist_median) + "\t" +

    str(undir_ref_1_inter_aa_dist_median) + "\t" +
    str(undir_ref_1_inter_ai_dist_median) + "\t" +
    str(undir_ref_1_inter_ii_dist_median) + "\t" +

    str(undir_ref_2_inter_aa_dist_median) + "\t" +
    str(undir_ref_2_inter_ai_dist_median) + "\t" +
    str(undir_ref_2_inter_ii_dist_median) +

    "\n"
)

tab_stream_stats_output.close()


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


### Create barplots for proportions of interactions within AA, AI and II
########################################################################

def create_barplot_for_enrichment_pair_tag_percentages(title, i_percentages, i_numbers):

    fig, ax = plt.subplots()

    xticklables_labels = ['DI', 'U', 'UR 1', 'UR 2']
    x = numpy.arange(len(xticklables_labels))  # the label locations
    width = 0.35  # the width of the bars
    ax.set_title(title)
    ax.set_xlabel('Interaction category')
    ax.set_ylabel('Percentages')
    ax.grid(zorder=0)
    ax.set_xticks(x)
    plt.ylim(0, max(i_percentages) + max(i_percentages)/10)
    ax.set_xticklabels(xticklables_labels)
    fig.tight_layout()

    rects = ax.bar(x, i_percentages, width, zorder=3)

    i = 0
    for rect in rects:
        height = rect.get_height()
        ax.annotate(i_numbers[i],
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        i += 1

    fig.tight_layout()

    return fig


pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name_barplots_interaction_enrichment_pair_tags)

aa_percentages = [round(100 * dir_inter_aa_num / dir_inter_num, 2), round(100 * undir_inter_aa_num / undir_inter_num, 2), round(100 * undir_ref_1_inter_aa_num / undir_ref_1_inter_num, 2), round(100 * undir_ref_2_inter_aa_num / undir_ref_2_inter_num, 2)]
aa_numbers = [dir_inter_aa_num, undir_inter_aa_num, undir_ref_1_inter_aa_num, undir_ref_2_inter_aa_num]

ai_percentages = [round(100 * dir_inter_ai_num / dir_inter_num, 2), round(100 * undir_inter_ai_num / undir_inter_num, 2), round(100 * undir_ref_1_inter_ai_num / undir_ref_1_inter_num, 2), round(100 * undir_ref_2_inter_ai_num / undir_ref_2_inter_num, 2)]
ai_numbers = [dir_inter_ai_num, undir_inter_ai_num, undir_ref_1_inter_ai_num, undir_ref_2_inter_ai_num]

ii_percentages = [round(100 * dir_inter_ii_num / dir_inter_num, 2), round(100 * undir_inter_ii_num / undir_inter_num, 2), round(100 * undir_ref_1_inter_ii_num / undir_ref_1_inter_num, 2), round(100 * undir_ref_2_inter_ii_num / undir_ref_2_inter_num, 2)]
ii_numbers = [dir_inter_ii_num, undir_inter_ii_num, undir_ref_1_inter_ii_num, undir_ref_2_inter_ii_num]

pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within AA", aa_percentages, aa_numbers))
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within AI", ai_percentages, ai_numbers))
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within II", ii_percentages, ii_numbers))


### Create barplots for proportions of interaction associated digests within AA, AI and II
##########################################################################################

a_percentages = [round(percentage_dir_a_dig, 2), round(percentage_undir_a_dig, 2), round(percentage_undir_ref_1_a_dig, 2), round(percentage_undir_ref_2_a_dig, 2)]
a_numbers = [dir_a_dig_num, undir_a_dig_num, undir_ref_1_a_dig_num, undir_ref_2_a_dig_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of active (A) digests involved in interaction categories", a_percentages, a_numbers))

pdf.close()