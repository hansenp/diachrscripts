#!/usr/bin/env python

"""
This script takes an enhanced interaction (EI) file (created with the script '05_define_di_uie_and_uii.py') with the
following interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions
      i. UIE - Exclusive undirected interactions. Subset of UI. No digest interacts with a digest involved in DI
      ii. UII - Inclusive undirected interactions. Subset of UI. UI without UIE

that are specified in the third column.

In addition, a BED file with the coordinates of the digests that were selected for enrichment must be specified. These
coordinates are used to assign one of the three enrichment states EE, EN and NN to each interaction, where 'E' stands
for 'enriched' and 'N' for 'not enriched'. In the EI file that will be generated, the enrichment status is written to
column 6, with the order of 'E' and 'N' indicating which of the two digests involved has been enriched.

Undirected reference interactions (UIR) that are comparable to DI with respect to enrichment status and read pair
numbers will be selected from UIE and UII.

The script implements a two pass approach:

   1. Iterate interactions and collect information about DI with respect to enrichment status of digests and read pair
   numbers.

   2. Iterate interactions a second time and use the collected information to select reference interactions that are
   comparable to DI.

The script implements the following approach for the selection of undirected reference interactions:

   Determine the exact distribution of read pair numbers within DI during the first pass and select the same number
   of undirected interactions for each read number during the second pass.

Directed interactions and undirected reference interactions are written to a file in enhanced interaction format with
the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_uir.tsv.gz'

Directed interactions and undirected interactions are written to a file in enhanced interaction format with
the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_ui.tsv.gz'

Column 3 of the created files contains the following interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions
   3. UIR - Undirected reference interactions

Finally, summary statistics about numbers of interactions within different enrichment state categories (EE, EN, NE and NN)
and corresponding read pair numbers are printed to the screen and a boxplot for the distribution of read pair numbers
will be created:

   '<OUT_PREFIX>_read_pair_number_boxplot.pdf'

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
parser.add_argument('--enriched-digests-file', help='BED file with digests selcted for target enrichment.')

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
enriched_digests_file = args.enriched_digests_file

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --enriched-digests-file: " + str(enriched_digests_file))

### Prepare variables
#####################

# Directed interactions
# ---------------------

# Dictionaries for read pair numbers of directed interactions within EE, EN and NN
di_ee_rp_dict = {}
di_en_rp_dict = {}
di_nn_rp_dict = {}

# Total number of interactions
di_num = 0

# Numbers of interactions within EE, EN and NN
di_ee_num = 0
di_en_num = 0
di_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
di_ee_rp_array = []
di_en_rp_array = []
di_nn_rp_array = []

# Sets for involved digests
di_e_dig_set = set()
di_n_dig_set = set()

# Undirected interactions
# -----------------------

# Total number of interactions
ui_num = 0

# Numbers of interactions within EE, EN and NN
ui_ee_num = 0
ui_en_num = 0
ui_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
ui_ee_rp_array = []
ui_en_rp_array = []
ui_nn_rp_array = []

# Sets for involved digests
ui_e_dig_set = set()
ui_n_dig_set = set()

# Undirected reference interactions
# ---------------------------------

# Total number of interactions
uir_num = 0

# Numbers of interactions within EE, EN and NN
uir_ee_num = 0
uir_en_num = 0
uir_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
uir_ee_rp_array = []
uir_en_rp_array = []
uir_nn_rp_array = []

# Sets for involved digests
uir_e_dig_set = set()
uir_n_dig_set = set()


### Read list with digests selected for enrichment
##################################################

print("[INFO] Reading list with digests selected for enrichment ...")
enriched_digests_set = set()
with open(enriched_digests_file, 'rt') as fp:
    line = fp.readline().rstrip()
    while line:

        # Parse line
        chr, sta, end = line.split('\t')

        # Add to set
        enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))

        # Go to next line
        line = fp.readline().rstrip()

fp.close()
print("\t[INFO] Read " + str(len(enriched_digests_set)) + " digests ...")
print("[INFO] ... done.")


### Prepare output files
########################

# Enhanced interaction file with interaction categories: DI and UIR
enhanced_interaction_file_output = out_prefix + "_enhanced_interaction_file_with_di_and_uir.tsv.gz"
enhanced_interaction_stream_output = gzip.open(enhanced_interaction_file_output, 'wt')

# PDF file with boxplots for distributions of read pair numbers
pdf_name_boxplots_read_pair_numbers = out_prefix + "_read_pair_number_boxplot.pdf"

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

        # Get digest coordinates from enhanced interaction file
        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        if coord_key_da in enriched_digests_set:
            enrichment_tag_dig_1 = 'E'
        else:
            enrichment_tag_dig_1 = 'N'

        if coord_key_db in enriched_digests_set:
            enrichment_tag_dig_2 = 'E'
        else:
            enrichment_tag_dig_2 = 'N'

        enrichment_pair_tag = enrichment_tag_dig_1 + enrichment_tag_dig_2

        # Collect read pair numbers for DI, UIE and UII
        if interaction_category == 'DI':

            di_num += 1

            if enrichment_pair_tag == 'EE':
                di_ee_num +=1
                di_ee_rp_array.append(rp_total)
                if rp_total not in di_ee_rp_dict:
                    di_ee_rp_dict[rp_total] = 1
                else:
                    di_ee_rp_dict[rp_total] +=1

            elif enrichment_pair_tag == 'EN' or  enrichment_pair_tag == 'NE':
                di_en_num += 1
                di_en_rp_array.append(rp_total)
                if rp_total not in di_en_rp_dict:
                    di_en_rp_dict[rp_total] = 1
                else:
                    di_en_rp_dict[rp_total] +=1

            elif enrichment_pair_tag == 'NN':
                di_nn_num += 1
                di_nn_rp_array.append(rp_total)
                if rp_total not in di_nn_rp_dict:
                    di_nn_rp_dict[rp_total] = 1
                else:
                    di_nn_rp_dict[rp_total] +=1

            # Collect enriched and not enriched digests
            if enrichment_tag_dig_1 == 'E':
                di_e_dig_set.add(coord_key_da)
            else:
                di_n_dig_set.add(coord_key_da)
            if enrichment_tag_dig_2 == 'E':
                di_e_dig_set.add(coord_key_db)
            else:
                di_n_dig_set.add(coord_key_db)

        elif interaction_category == 'UIE' or interaction_category == 'UII':

            ui_num += 1

            if enrichment_pair_tag == 'EE':
                ui_ee_num +=1
                ui_ee_rp_array.append(rp_total)

            elif enrichment_pair_tag == 'EN' or  enrichment_pair_tag == 'NE':
                ui_en_num += 1
                ui_en_rp_array.append(rp_total)

            elif enrichment_pair_tag == 'NN':
                ui_nn_num += 1
                ui_nn_rp_array.append(rp_total)

            # Collect enriched and not enriched digests
            if enrichment_tag_dig_1 == 'E':
                ui_e_dig_set.add(coord_key_da)
            else:
                ui_n_dig_set.add(coord_key_da)
            if enrichment_tag_dig_2 == 'E':
                ui_e_dig_set.add(coord_key_db)
            else:
                ui_n_dig_set.add(coord_key_db)

        else:
            print("[Error] Interaction category must be either \'DI\', \'UIE\' or \'UII\'!")
            exit(1)

        line = fp.readline()

fp.close()
print("\t[INFO] done ...")

if di_ee_num < 3 or di_en_num < 3 or di_nn_num < 3:

    print("[ERROR] Too few directed interactions within EE, EN or NN! Cannot select reference interactions.")

    print("\tdi_num\tdi_ee_num\tdi_en_num\tdi_nn_num")
    print("\t" + str(di_num) + "\t" + str(di_ee_num) + "\t" + str(di_en_num) + "\t" + str(di_nn_num))

    print("\tui_num\tui_ee_num\tui_en_num\tui_nn_num")
    print("\t" + str(ui_num) + "\t" + str(ui_ee_num) + "\t" + str(ui_en_num) + "\t" + str(ui_nn_num))

    tab_stream_warnings_output.close()
    tab_stream_stats_output.close()

    exit(0)


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

        # Get digest coordinates from enhanced interaction file
        coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
        coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)

        if coord_key_da in enriched_digests_set:
            enrichment_tag_dig_1 = 'E'
        else:
            enrichment_tag_dig_1 = 'N'

        if coord_key_db in enriched_digests_set:
            enrichment_tag_dig_2 = 'E'
        else:
            enrichment_tag_dig_2 = 'N'

        enrichment_pair_tag = enrichment_tag_dig_1 + enrichment_tag_dig_2

        # Write directed interaction to file
        if interaction_category == 'DI':
            enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, 'DI') + "\n")

        # Select undirected reference interactions from UII aand UIE using dictionaries for read pair numbers in DI
        if interaction_category == 'UIE' or interaction_category == 'UII':
            interaction_category_tag = 'UI'

            if enrichment_pair_tag == 'EE' and rp_total in di_ee_rp_dict and 0 < di_ee_rp_dict[rp_total]:
                di_ee_rp_dict[rp_total] -= 1
                uir_ee_rp_array.append(rp_total)
                uir_num +=1
                uir_ee_num += 1
                interaction_category_tag = 'UIR'

                uir_e_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                uir_e_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            elif (enrichment_pair_tag == 'EN' or enrichment_pair_tag == 'NE') and rp_total in di_en_rp_dict and 0 < di_en_rp_dict[rp_total]:
                di_en_rp_dict[rp_total] -= 1
                uir_en_rp_array.append(rp_total)
                uir_num += 1
                uir_en_num += 1
                interaction_category_tag = 'UIR'

                if enrichment_tag_dig_1 == 'E':
                    uir_e_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                else:
                    uir_n_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                if enrichment_tag_dig_2 == 'E':
                    uir_e_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))
                else:
                    uir_n_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            elif enrichment_pair_tag == 'NN' and rp_total in di_nn_rp_dict and 0 < di_nn_rp_dict[rp_total]:
                di_nn_rp_dict[rp_total] -= 1
                uir_nn_rp_array.append(rp_total)
                uir_num += 1
                uir_nn_num += 1
                interaction_category_tag = 'UIR'

                uir_n_dig_set.add(chr_a + "\t" + str(sta_a) + "\t" + str(end_a))
                uir_n_dig_set.add(chr_b + "\t" + str(sta_b) + "\t" + str(end_b))

            # Override interaction category tag in column 3 and write interaction to file
            enhanced_interaction_stream_output.write(diachrscripts_toolkit.set_interaction_category_in_enhanced_interaction_line(line, interaction_category_tag) + "\n")

        line = fp.readline()

enhanced_interaction_stream_output.close()
fp.close()
print("\t[INFO] done ...")

print("[INFO] Checking whether reference interactions are missing for some n")

tab_stream_warnings_output.write("EE_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_aa = 0
for x in di_ee_rp_dict:
    if 0 < di_ee_rp_dict[x]:
        n_missing_aa += di_ee_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_ee_rp_dict[x]) + "\n")

tab_stream_warnings_output.write("EN_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_ai = 0
for x in di_en_rp_dict:
    if 0 < di_en_rp_dict[x]:
        n_missing_ai += di_en_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_en_rp_dict[x]) + "\n")

tab_stream_warnings_output.write("NN_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_ii = 0
for x in di_nn_rp_dict:
    if 0 < di_nn_rp_dict[x]:
        n_missing_ii += di_nn_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_nn_rp_dict[x]) + "\n")

print("\t[INFO] Summary for EE, EN and NN")
print("\t\t[WARNING] For " + str(n_missing_aa) + " out of " + str(di_ee_num) + " directed EE interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_ai) + " out of " + str(di_en_num) + " directed EN interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_ii) + " out of " + str(di_nn_num) + " directed NN interactions no undirected reference interaction could be selected.")

tab_stream_warnings_output.close()
print("\t[INFO] done ...")


### Output some statistics
##########################

di_ee_rp_median = int(numpy.quantile(di_ee_rp_array, 0.50))
di_en_rp_median = int(numpy.quantile(di_en_rp_array, 0.50))
di_nn_rp_median = int(numpy.quantile(di_nn_rp_array, 0.50))

di_ee_percent = "{0:.2f}".format(100 * di_ee_num / di_num)
di_en_percent = "{0:.2f}".format(100 * di_en_num / di_num)
di_nn_percent = "{0:.2f}".format(100 * di_nn_num / di_num)

ui_ee_rp_median = int(numpy.quantile(ui_ee_rp_array, 0.50))
ui_en_rp_median = int(numpy.quantile(ui_en_rp_array, 0.50))
ui_nn_rp_median = int(numpy.quantile(ui_nn_rp_array, 0.50))

ui_ee_percent = "{0:.2f}".format(100 * ui_ee_num / ui_num)
ui_en_percent = "{0:.2f}".format(100 * ui_en_num / ui_num)
ui_nn_percent = "{0:.2f}".format(100 * ui_nn_num / ui_num)

uir_ee_rp_median = int(numpy.quantile(uir_ee_rp_array, 0.50))
uir_en_rp_median = int(numpy.quantile(uir_en_rp_array, 0.50))
uir_nn_rp_median = int(numpy.quantile(uir_nn_rp_array, 0.50))

uir_ee_percent = "{0:.2f}".format(100 * uir_ee_num / uir_num)
uir_en_percent = "{0:.2f}".format(100 * uir_en_num / uir_num)
uir_nn_percent = "{0:.2f}".format(100 * uir_nn_num / uir_num)

print()
print("Total number of directed interactions: " + str(di_num))
print("\tWithin 'EE': " + str(di_ee_num) + " (" + str(di_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_ee_rp_median))
print("\tWithin 'EN': " + str(di_en_num) + " (" + str(di_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_en_rp_median))
print("\tWithin 'NN': " + str(di_nn_num) + " (" + str(di_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_nn_rp_median))

print()
print("Number of exclusive and inclusive undirected interactions: " + str(ui_num))
print("\tWithin 'EE': " + str(ui_ee_num) + " (" + str(ui_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_ee_rp_median))
print("\tWithin 'EN': " + str(ui_en_num) + " (" + str(ui_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_en_rp_median))
print("\tWithin 'NN': " + str(ui_nn_num) + " (" + str(ui_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_nn_rp_median))

print()
print("Number of exclusive and inclusive undirected reference (2) interactions: " + str(uir_num))
print("\tWithin 'EE': " + str(uir_ee_num) + " (" + str(uir_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_ee_rp_median))
print("\tWithin 'EN': " + str(uir_en_num) + " (" + str(uir_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_en_rp_median))
print("\tWithin 'NN': " + str(uir_nn_num) + " (" + str(uir_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_nn_rp_median))

dir_a_dig_num = len(di_e_dig_set)                   # Number of active digests involved in directed interactions
dir_i_dig_num = len(di_n_dig_set)                   # Number of inactive digests involved in directed interactions
percentage_dir_a_dig = 100 * dir_a_dig_num / (dir_a_dig_num + dir_i_dig_num)

undir_a_dig_num = len(ui_e_dig_set)               # Number of active digests involved in undirected interactions
undir_i_dig_num = len(ui_n_dig_set)               # Number of inactive digests involved in undirected interactions
percentage_undir_a_dig = 100 * undir_a_dig_num / (undir_a_dig_num + undir_i_dig_num)

undir_ref_a_dig_num = len(uir_e_dig_set)   # Number of active digests involved in undirected interactions
undir_ref_i_dig_num = len(uir_n_dig_set)   # Number of inactive digests involved in undirected interactions
percentage_undir_ref_a_dig = 100 * undir_ref_a_dig_num / (undir_ref_a_dig_num + undir_ref_i_dig_num)

print()
print("Number of active (E) digests involved in DI: " + str(dir_a_dig_num) + " (" + "{0:.2f}".format(percentage_dir_a_dig) + "%)")
print("Number of inactive (N) digests involved in DI: " + str(dir_i_dig_num))
print()
print("Number of active (E) digests involved in U: " + str(undir_a_dig_num) + " (" + "{0:.2f}".format(percentage_undir_a_dig) + "%)")
print("Number of inactive (N) digests involved in U: " + str(undir_i_dig_num))
print()
print("Number of active (E) digests involved in UR2: " + str(undir_ref_a_dig_num) + " (" + "{0:.2f}".format(percentage_undir_ref_a_dig) + "%)")
print("Number of inactive (N) digests involved in UR2: " + str(undir_ref_i_dig_num))
print()

print()

tab_stream_stats_output.write(

    "out_prefix" + "\t" +                               # Prefix for output

    "di_num" + "\t" +                            # Total number of directed interactions

    "di_ee_num" + "\t" +                         # Number of directed interactions within AA
    "di_en_num" + "\t" +                         # Number of directed interactions within AI
    "di_nn_num" + "\t" +                         # Number of directed interactions within II

    "di_ee_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within AA
    "di_en_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within AI
    "di_nn_rp_median" + "\t" +                   # Median number of read pairs in directed interactions within II

    "ui_num" + "\t" +                          # Total number of undirected interactions

    "ui_ee_num" + "\t" +                       # Number of undirected interactions within AA
    "ui_en_num" + "\t" +                       # Number of undirected interactions within AI
    "ui_nn_num" + "\t" +                       # Number of undirected interactions within II

    "ui_ee_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within AA
    "ui_en_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within AI
    "ui_nn_rp_median" + "\t" +                 # Median number of read pairs in undirected interactions within II

    "uir_num" + "\t" +                    # Total number of undirected interactions (exact)

    "uir_ee_num" + "\t" +                 # Number of undirected reference interactions within AA
    "uir_en_num" + "\t" +                 # Number of undirected reference interactions within AI
    "uir_nn_num" + "\t" +                 # Number of undirected reference interactions within II

    "n_missing_aa" + "\t" +                             # Number of missing reference interactions within AA
    "n_missing_ai" + "\t" +                             # Number of missing reference interactions within AI
    "n_missing_ii" + "\t" +                             # Number of missing reference interactions within II

    "uir_ee_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AA
    "uir_en_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within AI
    "uir_nn_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within II

    "dir_a_dig_num" + "\t" +                            # Number of active digests involved in directed interactions
    "dir_i_dig_num" + "\t" +                            # Number of inactive digests involved in directed interactions

    "undir_a_dig_num" + "\t" +                          # Number of active digests involved in undirected interactions
    "undir_i_dig_num" + "\t" +                          # Number of inactive digests involved in undirected interactions

    "undir_ref_a_dig_num" + "\t" +                    # Number of active digests involved in undirected reference interactions (exact)
    "undir_ref_i_dig_num" +                           # Number of inactive digests involved in undirected reference interactions (exact)

    "\n"
)
tab_stream_stats_output.write(

    str(out_prefix) + "\t" +

    str(di_num) + "\t" +

    str(di_ee_num) + "\t" +
    str(di_en_num) + "\t" +
    str(di_nn_num) + "\t" +

    str(di_ee_rp_median) + "\t" +
    str(di_en_rp_median) + "\t" +
    str(di_nn_rp_median) + "\t" +

    str(ui_num) + "\t" +

    str(ui_ee_num) + "\t" +
    str(ui_en_num) + "\t" +
    str(ui_nn_num) + "\t" +

    str(ui_ee_rp_median) + "\t" +
    str(ui_en_rp_median) + "\t" +
    str(ui_nn_rp_median) + "\t" +

    str(uir_num) + "\t" +

    str(uir_ee_num) + "\t" +
    str(uir_en_num) + "\t" +
    str(uir_nn_num) + "\t" +

    str(n_missing_aa) + "\t" +
    str(n_missing_ai) + "\t" +
    str(n_missing_ii) + "\t" +

    str(uir_ee_rp_median) + "\t" +
    str(uir_en_rp_median) + "\t" +
    str(uir_nn_rp_median) + "\t" +

    str(dir_a_dig_num) + "\t" +
    str(dir_i_dig_num) + "\t" +

    str(undir_a_dig_num) + "\t" +
    str(undir_i_dig_num) + "\t" +

    str(undir_ref_a_dig_num) + "\t" +
    str(undir_ref_i_dig_num) +

    "\n"
)

tab_stream_stats_output.close()


### Create boxplots for read pair numbers
#########################################

plt.rcParams.update({'font.size':15})

data = [
    di_ee_rp_array, di_en_rp_array, di_nn_rp_array,
    ui_ee_rp_array,ui_en_rp_array, ui_nn_rp_array,
    uir_ee_rp_array, uir_en_rp_array, uir_nn_rp_array
]
labels = ['EE', 'EN', 'NN', 'EE', 'EN', 'NN','EE', 'EN', 'NN']
fig1, ax1 = plt.subplots()
ax1.set_title('Distributions of read pair numbers for DI, U, and UR2 (' + out_prefix + ')')
ax1.boxplot(data, showfliers=False, labels=labels)

box = ax1.boxplot(data, showfliers=False, labels=labels, patch_artist=True)

colors = ['orange', 'orange', 'orange', 'darkgray', 'darkgray', 'darkgray', 'lightblue', 'lightblue', 'lightblue']

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

lab = "DI - AA: " + str(di_ee_num) + " (" + di_ee_percent + "%), AI: " + str(di_en_num) + " (" + di_en_percent + "%), II: " + str(di_nn_num) + " (" + di_nn_percent + "%)"
di_patch = mpatches.Patch(color='orange', label=lab)
lab = "U - AA: " + str(ui_ee_num) + " (" + ui_ee_percent + "%), AI: " + str(ui_en_num) + " (" + ui_en_percent + "%), II: " + str(ui_nn_num) + " (" + ui_nn_percent + "%)"
u_patch = mpatches.Patch(color='darkgray', label=lab)
lab = "UR 2 - AA: " + str(uir_ee_num) + " (" + uir_ee_percent + "%), AI: " + str(uir_en_num) + " (" + uir_en_percent + "%), II: " + str(uir_nn_num) + " (" + uir_nn_percent + "%)"
ur2_patch = mpatches.Patch(color='lightblue', label=lab)
plt.legend(handles=[di_patch, u_patch, ur2_patch])
ax1.set_xlabel('Enrichment state pair tag')
ax1.set_ylabel('Read pair number')
plt.grid(True)
fig1.set_size_inches(7,5)
plt.savefig(pdf_name_boxplots_read_pair_numbers, bbox_inches='tight')
plt.close()


### Create barplots for proportions of interactions within EE, EN and NN
########################################################################

def create_barplot_for_enrichment_pair_tag_percentages(title, i_percentages, i_numbers):

    fig, ax = plt.subplots()

    xticklables_labels = ['DI', 'U', 'UR']
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

aa_percentages = [round(100 * di_ee_num / di_num, 2), round(100 * ui_ee_num / ui_num, 2), round(100 * uir_ee_num / uir_num, 2)]
aa_numbers = [di_ee_num, ui_ee_num, uir_ee_num]

ai_percentages = [round(100 * di_en_num / di_num, 2), round(100 * ui_en_num / ui_num, 2), round(100 * uir_en_num / uir_num, 2)]
ai_numbers = [di_en_num, ui_en_num, uir_en_num]

ii_percentages = [round(100 * di_nn_num / di_num, 2), round(100 * ui_nn_num / ui_num, 2), round(100 * uir_nn_num / uir_num, 2)]
ii_numbers = [di_nn_num, ui_nn_num, uir_nn_num]

pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within EE", aa_percentages, aa_numbers))
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within EN", ai_percentages, ai_numbers))
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within NN", ii_percentages, ii_numbers))


### Create barplots for proportions of interaction associated digests within AA, AI and II
##########################################################################################

a_percentages = [round(percentage_dir_a_dig, 2), round(percentage_undir_a_dig, 2), round(percentage_undir_ref_a_dig, 2)]
a_numbers = [dir_a_dig_num, undir_a_dig_num, undir_ref_a_dig_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of active (E) digests involved in interaction categories", a_percentages, a_numbers))

pdf.close()