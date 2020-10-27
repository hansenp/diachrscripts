#!/usr/bin/env python

"""
This script takes an enhanced interaction (EI) file (created with the script '05_define_di_uie_and_uii.py') with the
following interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions
      i. UIE - Exclusive undirected interactions. Subset of UI. No digest interacts with a digest involved in DI
      ii. UII - Inclusive undirected interactions. Subset of UI. UI without UIE

that are specified in the third column of the EI file.

In addition, a BED file with the coordinates of the digests that were selected for enrichment must be specified. These
coordinates are used to assign one of the four enrichment states EE, EN, NE and NN to each interaction, where 'E' stands
for 'enriched' and 'N' for 'not enriched'. In the EI file that will be generated, the enrichment status is written to
column 6, with the order of 'E' and 'N' indicating which of the two digests involved has been enriched.

Undirected reference interactions (UIR) that are comparable to DI with respect to the distribution of read pair numbers
per interaction within the enrichment states EE, EN and NN will be selected from UIE and UII.

The script implements a two pass approach:

   1. Iterate interactions and determine the distribution of read pair numbers per interaction within DI.

   2. Iterate interactions a second time and select a set of undirected reference interactions that are comparable to
    DI with respect to the distribution of read pair numbers per interaction.

Directed interactions and undirected reference interactions are written to a file in enhanced interaction format with
the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_uir.tsv.gz'

Directed interactions and all undirected interactions including reference interactions are written to a second file with
the name:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz'

Column 3 of the created files contains the following tags for the interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions
   3. UIR - Undirected reference interactions

Column 6 of the created files contains the following the following tags for the enrichment states of interactions:

   1. EE - Both digests enriched
   2. NE - Second digests enriched
   3. EN - First digests enriched
   4. NN - No digest enriched

Finally, summary statistics are written to a tab separated file:

   '<OUT_PREFIX>_stats.tsv'

and partly also to the screen. For some directed interactions, there is no undirected interaction with the same number
of read pairs (typically for large read pair numbers) and, therefore, no reference interaction can be selected.
Statistics about this will be reported in the file:

   '<OUT_PREFIX>_warnings.tsv'

Furthermore, a boxplot with the distribution of read pair numbers will be created:

   '<OUT_PREFIX>_read_pair_number_boxplot.pdf'

Finally, barplots for the proportions of interactions within DI, UI and UIR:

   '<OUT_PREFIX>_interaction_enrichment_pair_tags_barplot.pdf'

Two EI files for downstream analyzes will be created, one file containing DI and UIR interactions only:

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_and_uir.tsv.gz'

and another EI file containing all input interactions, but with overwritten tags for interaction categories (column 3)
and enrichment states (column 6):

   '<OUT_PREFIX>_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz

Explain option '--respect-left-right'

"""


import argparse
import gzip
import diachrscripts_toolkit
from diachr2 import EnhancedInteraction, EnhancedInteractionParser
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
parser.add_argument('--respect-left-right', help='When choosing the reference interactions, treat NE and EN as separate categories.', action='store_true', default=False)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
enriched_digests_file = args.enriched_digests_file
respect_left_right = args.respect_left_right

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)
print("\t[INFO] --enriched-digests-file: " + str(enriched_digests_file))
print("\t[INFO] --respect-left-right: " + str(respect_left_right))


### Prepare variables
#####################

# Directed interactions
# ---------------------

# Dictionaries for read pair numbers of directed interactions within EE, EN and NN
di_ee_rp_dict = {}
di_ne_rp_dict = {} #
di_en_rp_dict = {}
di_neen_rp_dict = {}
di_nn_rp_dict = {}

# Total number of interactions
di_num = 0

# Numbers of interactions within EE, EN and NN
di_ee_num = 0
di_en_num = 0
di_ne_num = 0
di_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
di_ee_rp_array = []
di_en_rp_array = []
di_ne_rp_array = []
di_nn_rp_array = []

# Undirected interactions
# -----------------------

# Total number of interactions
ui_num = 0

# Numbers of interactions within EE, EN and NN
ui_ee_num = 0
ui_en_num = 0
ui_ne_num = 0
ui_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
ui_ee_rp_array = []
ui_en_rp_array = []
ui_ne_rp_array = []
ui_nn_rp_array = []

# Undirected reference interactions
# ---------------------------------

# Total number of interactions
uir_num = 0

# Numbers of interactions within EE, EN and NN
uir_ee_num = 0
uir_en_num = 0
uir_ne_num = 0
uir_nn_num = 0

# Arrays for read pair numbers within EE, EN and NN
uir_ee_rp_array = []
uir_en_rp_array = []
uir_ne_rp_array = []
uir_nn_rp_array = []


### Read list with digests selected for enrichment
##################################################

print("[INFO] Reading list with digests selected for enrichment ...")
enriched_digests_set = set()
with open(enriched_digests_file, 'rt') as fp:
    for line in fp:
        chr, sta, end = line.rstrip().split('\t')
        enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))

print("\t[INFO] Read " + str(len(enriched_digests_set)) + " digests ...")
print("[INFO] ... done.")


### Prepare output files
########################

# EI file with interaction categories: DI and UIR only
di_uir_enhanced_interaction_stream_output = gzip.open(out_prefix + "_enhanced_interaction_file_with_di_and_uir.tsv.gz", 'wt')

# EI file with all input interactions, but with overwritten tags for interaction categories (column 3) and enrichment states (column 6)
di_ui_uir_enhanced_interaction_stream_output = gzip.open(out_prefix + "_enhanced_interaction_file_with_di_ui_and_uir.tsv.gz", 'wt')

# Prepare stream for output of statistics
tab_stream_stats_output = open(out_prefix + "_stats.tsv", 'wt')

# Prepare stream for warnings
tab_stream_warnings_output = open(out_prefix + "_warnings.tsv", 'wt')

# PDF file with boxplots for distributions of read pair numbers
pdf_name_boxplots_read_pair_numbers = out_prefix + "_read_pair_number_boxplot.pdf"

# PDF file with three barplots showing the proportions of interactions within the enrichment pair categories EE, EN, NE and NN
pdf_name_barplots_interaction_enrichment_pair_tags = out_prefix + "_interaction_enrichment_pair_tags_barplot.pdf"


### 1st pass: Determine distribution of read pairs per interaction for DI
#########################################################################

print("[INFO] 1st pass: Collect information about DI ...")

print("\t[INFO] Reading enhanced interaction file ...")
#with gzip.open(enhanced_interaction_file, 'rt') as fp:
ie_parser = EnhancedInteractionParser(input_path)
ei_list = ie_parser.parse()
print("[INFO] Extracted %d interaction lines" % len(ei_list))

n_interaction_total = 0
#    for line in fp:
for ei in ei_list:

    # Report progress
    n_interaction_total += 1
    if n_interaction_total % 1000000 == 0:
        print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

    # Parse enhanced interactions line
    #chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
    #    diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)
    enrichment_pair_tag = ei.enrichment_pair_tag
    interaction_category = ei.interaction_category
    rp_total = ei.rp_total

    # Get digest coordinates from enhanced interaction file
    #coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
    #coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)
    coord_key_da = ei.coordinate_key_a
    coord_key_db = ei.coordinate_key_b

    # Get enrichment status of digests from file for enriched digests
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

        elif enrichment_pair_tag == 'EN':
            di_en_num += 1
            di_en_rp_array.append(rp_total)
            if rp_total not in di_en_rp_dict:
                di_en_rp_dict[rp_total] = 1
            else:
                di_en_rp_dict[rp_total] +=1
            if rp_total not in di_neen_rp_dict:
                di_neen_rp_dict[rp_total] = 1
            else:
                di_neen_rp_dict[rp_total] += 1

        elif enrichment_pair_tag == 'NE':
            di_ne_num += 1
            di_ne_rp_array.append(rp_total)
            if rp_total not in di_ne_rp_dict: #
                di_ne_rp_dict[rp_total] = 1 #
            else:
                di_ne_rp_dict[rp_total] +=1 #
            if rp_total not in di_neen_rp_dict:
                di_neen_rp_dict[rp_total] = 1
            else:
                di_neen_rp_dict[rp_total] += 1

        elif enrichment_pair_tag == 'NN':
            di_nn_num += 1
            di_nn_rp_array.append(rp_total)
            if rp_total not in di_nn_rp_dict:
                di_nn_rp_dict[rp_total] = 1
            else:
                di_nn_rp_dict[rp_total] +=1

    # Count all undirected interactions and collect read pair numbers
    elif interaction_category == 'UIE' or interaction_category == 'UII':

        ui_num += 1

        if enrichment_pair_tag == 'EE':
            ui_ee_num +=1
            ui_ee_rp_array.append(rp_total)

        elif enrichment_pair_tag == 'EN':
            ui_en_num += 1
            ui_en_rp_array.append(rp_total)

        elif enrichment_pair_tag == 'NE':
            ui_ne_num += 1
            ui_ne_rp_array.append(rp_total)

        elif enrichment_pair_tag == 'NN':
            ui_nn_num += 1
            ui_nn_rp_array.append(rp_total)

    else:
        print("[Error] Interaction category must be either \'DI\', \'UIE\' or \'UII\'!")
        exit(1)

print("\t[INFO] done ...")

if di_ee_num < 3 or di_en_num < 3 or di_nn_num < 3:

    print("[ERROR] Too few directed interactions within EE, EN or NN! Cannot select reference interactions.")

    print("\tdi_num\tdi_ee_num\tdi_en_num\tdi_ne_num\tdi_nn_num")
    print("\t" + str(di_num) + "\t" + str(di_ee_num) + "\t" + str(di_en_num) + "\t" + str(di_ne_num) + "\t" + str(di_nn_num))

    print("\tui_num\tui_ee_num\tui_en_num\tui_ne_num\tui_nn_num")
    print("\t" + str(ui_num) + "\t" + str(ui_ee_num) + "\t" + str(ui_en_num) + "\t" + str(ui_ne_num) + "\t" + str(ui_nn_num))

    exit(0)


### 2nd pass: Select undirected reference interactions
######################################################
print("[INFO] 2nd pass: Select undirected reference interactions ...")

print("\t[INFO] Iterating enhanced interaction file ...")
#with gzip.open(enhanced_interaction_file, 'rt') as fp:
ie_parser = EnhancedInteractionParser(input_path)
ei_list = ie_parser.parse()
print("[INFO] 2nd pass Extracted %d interaction lines" % len(ei_list))

n_interaction_total = 0
#for line in fp:
for ei in ei_list:
    n_interaction_total += 1
    if n_interaction_total % 1000000 == 0:
        print("\t\t[INFO]", n_interaction_total, "interactions processed ...")

    # Parse enhanced interactions line
    # chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
    #        diachrscripts_toolkit.parse_enhanced_interaction_line_with_gene_symbols(line)
    enrichment_pair_tag = ei.enrichment_pair_tag
    interaction_category = ei.interaction_category
    rp_total = ei.rp_total
    # Get digest coordinates from enhanced interaction file
    # coord_key_da = chr_a + '\t' + str(sta_a) + '\t' + str(end_a)
    # coord_key_db = chr_b + '\t' + str(sta_b) + '\t' + str(end_b)
    coord_key_da = ei.coordinate_key_a
    coord_key_db = ei.coordinate_key_b

    # Get enrichment status of digests from file for enriched digests
    if coord_key_da in enriched_digests_set:
        enrichment_tag_dig_1 = 'E'
    else:
        enrichment_tag_dig_1 = 'N'

    if coord_key_db in enriched_digests_set:
        enrichment_tag_dig_2 = 'E'
    else:
        enrichment_tag_dig_2 = 'N'

    enrichment_pair_tag = enrichment_tag_dig_1 + enrichment_tag_dig_2

    # Set tag for interaction category
    if interaction_category == 'DI':
        interaction_category_tag = 'DI'

    # Select undirected reference interactions from UII and UIE using dictionaries for read pair numbers in DI
    if interaction_category == 'UIE' or interaction_category == 'UII':

        interaction_category_tag = 'UI'

        if enrichment_pair_tag == 'EE' and rp_total in di_ee_rp_dict and 0 < di_ee_rp_dict[rp_total]:
            di_ee_rp_dict[rp_total] -= 1
            uir_ee_rp_array.append(rp_total)
            uir_num +=1
            uir_ee_num += 1
            interaction_category_tag = 'UIR'

        elif (enrichment_pair_tag == 'EN') or (enrichment_pair_tag == 'NE'):

            if respect_left_right:

                if (enrichment_pair_tag == 'EN') and rp_total in di_en_rp_dict and 0 < di_en_rp_dict[rp_total]:
                    di_en_rp_dict[rp_total] -= 1
                    uir_en_rp_array.append(rp_total)
                    uir_num += 1
                    uir_en_num += 1
                    interaction_category_tag = 'UIR'

                elif (enrichment_pair_tag == 'NE') and rp_total in di_ne_rp_dict and 0 < di_ne_rp_dict[rp_total]:
                    di_ne_rp_dict[rp_total] -= 1
                    uir_ne_rp_array.append(rp_total)
                    uir_num += 1
                    uir_ne_num += 1
                    interaction_category_tag = 'UIR'

            else:

                if (enrichment_pair_tag == 'EN'):

                    if rp_total in di_neen_rp_dict and 0 < di_neen_rp_dict[rp_total]:
                        di_neen_rp_dict[rp_total] -= 1
                        uir_en_rp_array.append(rp_total)
                        uir_num += 1
                        uir_en_num += 1
                        interaction_category_tag = 'UIR'

                elif (enrichment_pair_tag == 'NE'):

                    if rp_total in di_neen_rp_dict and 0 < di_neen_rp_dict[rp_total]:
                        di_neen_rp_dict[rp_total] -= 1
                        uir_ne_rp_array.append(rp_total)
                        uir_num += 1
                        uir_ne_num += 1
                        interaction_category_tag = 'UIR'

        elif enrichment_pair_tag == 'NN' and rp_total in di_nn_rp_dict and 0 < di_nn_rp_dict[rp_total]:
            di_nn_rp_dict[rp_total] -= 1
            uir_nn_rp_array.append(rp_total)
            uir_num += 1
            uir_nn_num += 1
            interaction_category_tag = 'UIR'

    # Override interaction category tag in column 3 and write interaction to files
    #line = diachrscripts_toolkit.set_column_in_enhanced_interaction_line(line, 6, enrichment_pair_tag)
    #line = diachrscripts_toolkit.set_column_in_enhanced_interaction_line(line, 3, interaction_category_tag)
    ei.set_enrichment_pair_tag(enrichment_pair_tag)
    ei.set_interaction_category(interaction_category_tag)
    di_ui_uir_enhanced_interaction_stream_output.write(line + "\n")
    if interaction_category_tag == 'DI' or interaction_category_tag == 'UIR':
        di_uir_enhanced_interaction_stream_output.write(line + "\n")


di_uir_enhanced_interaction_stream_output.close()
di_ui_uir_enhanced_interaction_stream_output.close()

print("\t[INFO] done ...")

print("[INFO] Checking whether reference interactions are missing for some n")

tab_stream_warnings_output.write("EE_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_ee = 0
for x in di_ee_rp_dict:
    if 0 < di_ee_rp_dict[x]:
        n_missing_ee += di_ee_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_ee_rp_dict[x]) + "\n")

tab_stream_warnings_output.write("EN_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_en = 0
for x in di_en_rp_dict:
    if 0 < di_en_rp_dict[x]:
        n_missing_en += di_en_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_en_rp_dict[x]) + "\n")

tab_stream_warnings_output.write("NN_INTERACTIONS" + "\n")
tab_stream_warnings_output.write("READ_PAIRS" + "\t" + "MISSING_REFERENCE_INTERACTIONS"  + "\n")
n_missing_nn = 0
for x in di_nn_rp_dict:
    if 0 < di_nn_rp_dict[x]:
        n_missing_nn += di_nn_rp_dict[x]
        tab_stream_warnings_output.write(str(x) + "\t" + str(di_nn_rp_dict[x]) + "\n")

print("\t[INFO] Summary for EE, EN and NN")
print("\t\t[WARNING] For " + str(n_missing_ee) + " out of " + str(di_ee_num) + " directed EE interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_en) + " out of " + str(di_en_num) + " directed EN or NE interactions no undirected reference interaction could be selected.")
print("\t\t[WARNING] For " + str(n_missing_nn) + " out of " + str(di_nn_num) + " directed NN interactions no undirected reference interaction could be selected.")

tab_stream_warnings_output.close()
print("\t[INFO] done ...")


### Output some statistics
##########################

di_ee_rp_median = int(numpy.quantile(di_ee_rp_array, 0.50))
di_en_rp_median = int(numpy.quantile(di_en_rp_array, 0.50))
di_ne_rp_median = int(numpy.quantile(di_ne_rp_array, 0.50))
di_nn_rp_median = int(numpy.quantile(di_nn_rp_array, 0.50))

di_ee_percent = "{0:.2f}".format(100 * di_ee_num / di_num)
di_en_percent = "{0:.2f}".format(100 * di_en_num / di_num)
di_ne_percent = "{0:.2f}".format(100 * di_ne_num / di_num)
di_nn_percent = "{0:.2f}".format(100 * di_nn_num / di_num)

ui_ee_rp_median = int(numpy.quantile(ui_ee_rp_array, 0.50))
ui_en_rp_median = int(numpy.quantile(ui_en_rp_array, 0.50))
ui_ne_rp_median = int(numpy.quantile(ui_ne_rp_array, 0.50))
ui_nn_rp_median = int(numpy.quantile(ui_nn_rp_array, 0.50))

ui_ee_percent = "{0:.2f}".format(100 * ui_ee_num / ui_num)
ui_en_percent = "{0:.2f}".format(100 * ui_en_num / ui_num)
ui_ne_percent = "{0:.2f}".format(100 * ui_ne_num / ui_num)
ui_nn_percent = "{0:.2f}".format(100 * ui_nn_num / ui_num)

uir_ee_rp_median = int(numpy.quantile(uir_ee_rp_array, 0.50))
uir_en_rp_median = int(numpy.quantile(uir_en_rp_array, 0.50))
uir_ne_rp_median = int(numpy.quantile(uir_ne_rp_array, 0.50))
uir_nn_rp_median = int(numpy.quantile(uir_nn_rp_array, 0.50))

uir_ee_percent = "{0:.2f}".format(100 * uir_ee_num / uir_num)
uir_en_percent = "{0:.2f}".format(100 * uir_en_num / uir_num)
uir_ne_percent = "{0:.2f}".format(100 * uir_ne_num / uir_num)
uir_nn_percent = "{0:.2f}".format(100 * uir_nn_num / uir_num)

print()
print("[INFO] Total number of directed interactions: " + str(di_num))

print("\tWithin 'EE': " + str(di_ee_num) + " (" + str(di_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_ee_rp_median))

print("\tWithin 'EN': " + str(di_en_num) + " (" + str(di_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_en_rp_median))

print("\tWithin 'NE': " + str(di_ne_num) + " (" + str(di_ne_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_ne_rp_median))

print("\tWithin 'NN': " + str(di_nn_num) + " (" + str(di_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(di_nn_rp_median))

print()
print("[INFO] Number of undirected interactions: " + str(ui_num))
print("\tWithin 'EE': " + str(ui_ee_num) + " (" + str(ui_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_ee_rp_median))

print("\tWithin 'EN': " + str(ui_en_num) + " (" + str(ui_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_en_rp_median))

print("\tWithin 'NE': " + str(ui_ne_num) + " (" + str(ui_ne_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_ne_rp_median))

print("\tWithin 'NN': " + str(ui_nn_num) + " (" + str(ui_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(ui_nn_rp_median))

print()
print("[INFO] Number of undirected reference interactions: " + str(uir_num))
print("\tWithin 'EE': " + str(uir_ee_num) + " (" + str(uir_ee_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_ee_rp_median))

print("\tWithin 'EN': " + str(uir_en_num) + " (" + str(uir_en_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_en_rp_median))

print("\tWithin 'NE': " + str(uir_ne_num) + " (" + str(uir_ne_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_ne_rp_median))

print("\tWithin 'NN': " + str(uir_nn_num) + " (" + str(uir_nn_percent) + "%)")
print("\t\tMedian number of read pairs: " + str(uir_nn_rp_median))

tab_stream_stats_output.write(

    "out_prefix" + "\t" +                 # Prefix for output

    "di_num" + "\t" +                     # Total number of directed interactions

    "di_ee_num" + "\t" +                  # Number of directed interactions within EE
    "di_en_num" + "\t" +                  # Number of directed interactions within EN
    "di_ne_num" + "\t" +                  # Number of directed interactions within NE
    "di_nn_num" + "\t" +                  # Number of directed interactions within NN

    "di_ee_rp_median" + "\t" +            # Median number of read pairs in directed interactions within EE
    "di_en_rp_median" + "\t" +            # Median number of read pairs in directed interactions within EN
    "di_ne_rp_median" + "\t" +            # Median number of read pairs in directed interactions within NE
    "di_nn_rp_median" + "\t" +            # Median number of read pairs in directed interactions within NN

    "ui_num" + "\t" +                     # Total number of undirected interactions

    "ui_ee_num" + "\t" +                  # Number of undirected interactions within EE
    "ui_en_num" + "\t" +                  # Number of undirected interactions within EN
    "ui_ne_num" + "\t" +                  # Number of undirected interactions within NE
    "ui_nn_num" + "\t" +                  # Number of undirected interactions within NN

    "ui_ee_rp_median" + "\t" +            # Median number of read pairs in undirected interactions within EE
    "ui_en_rp_median" + "\t" +            # Median number of read pairs in undirected interactions within EN
    "ui_ne_rp_median" + "\t" +            # Median number of read pairs in undirected interactions within NE
    "ui_nn_rp_median" + "\t" +            # Median number of read pairs in undirected interactions within NN

    "uir_num" + "\t" +                    # Total number of undirected interactions (exact)

    "uir_ee_num" + "\t" +                 # Number of undirected reference interactions within EE
    "uir_en_num" + "\t" +                 # Number of undirected reference interactions within EN
    "uir_ne_num" + "\t" +                 # Number of undirected reference interactions within NE
    "uir_nn_num" + "\t" +                 # Number of undirected reference interactions within NN

    "n_missing_ee" + "\t" +               # Number of missing reference interactions within EE
    "n_missing_en" + "\t" +               # Number of missing reference interactions within EN or NE
    "n_missing_nn" + "\t" +               # Number of missing reference interactions within NN

    "uir_ee_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within EE
    "uir_en_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within EN
    "uir_ne_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within NE
    "uir_nn_rp_median" + "\t" +           # Median number of read pairs in undirected reference interactions within NN

    "\n"
)
tab_stream_stats_output.write(

    str(out_prefix) + "\t" +

    str(di_num) + "\t" +

    str(di_ee_num) + "\t" +
    str(di_en_num) + "\t" +
    str(di_ne_num) + "\t" +
    str(di_nn_num) + "\t" +

    str(di_ee_rp_median) + "\t" +
    str(di_en_rp_median) + "\t" +
    str(di_ne_rp_median) + "\t" +
    str(di_nn_rp_median) + "\t" +

    str(ui_num) + "\t" +

    str(ui_ee_num) + "\t" +
    str(ui_en_num) + "\t" +
    str(ui_ne_num) + "\t" +
    str(ui_nn_num) + "\t" +

    str(ui_ee_rp_median) + "\t" +
    str(ui_en_rp_median) + "\t" +
    str(ui_ne_rp_median) + "\t" +
    str(ui_nn_rp_median) + "\t" +

    str(uir_num) + "\t" +

    str(uir_ee_num) + "\t" +
    str(uir_en_num) + "\t" +
    str(uir_ne_num) + "\t" +
    str(uir_nn_num) + "\t" +

    str(n_missing_ee) + "\t" +
    str(n_missing_en) + "\t" +
    str(n_missing_nn) + "\t" +

    str(uir_ee_rp_median) + "\t" +
    str(uir_en_rp_median) + "\t" +
    str(uir_ne_rp_median) + "\t" +
    str(uir_nn_rp_median) +

    "\n"
)
tab_stream_stats_output.close()


### Create boxplots for read pair numbers
#########################################

plt.rcParams.update({'font.size':15})

data = [
    di_ee_rp_array, di_en_rp_array, di_ne_rp_array, di_nn_rp_array,
    ui_ee_rp_array, ui_en_rp_array, ui_ne_rp_array, ui_nn_rp_array,
    uir_ee_rp_array, uir_en_rp_array, uir_ne_rp_array, uir_nn_rp_array
]
labels = ['EE', 'EN', 'NE', 'NN', 'EE', 'EN', 'NE', 'NN','EE', 'EN', 'NE', 'NN']
fig1, ax1 = plt.subplots()
ax1.set_title('Distributions of read pair numbers for DI, UI, and UIR (' + out_prefix + ')', fontsize=10)
ax1.boxplot(data, showfliers=False, labels=labels)

box = ax1.boxplot(data, showfliers=False, labels=labels, patch_artist=True)

colors = ['orange', 'orange', 'orange', 'orange', 'darkgray', 'darkgray', 'darkgray', 'darkgray', 'lightblue', 'lightblue', 'lightblue', 'lightblue']

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

lab = "DI - EE: " + str(di_ee_num) + " (" + di_ee_percent + "%), EN: " + str(di_en_num) + " (" + di_en_percent + "%), NE: " + str(di_ne_num) + " (" + di_ne_percent + "%), NN: " + str(di_nn_num) + " (" + di_nn_percent + "%)"
di_patch = mpatches.Patch(color='orange', label=lab)
lab = "UI - EE: " + str(ui_ee_num) + " (" + ui_ee_percent + "%), EN: " + str(ui_en_num) + " (" + ui_en_percent + "%), NE: " + str(ui_ne_num) + " (" + ui_ne_percent + "%), NN: " + str(ui_nn_num) + " (" + ui_nn_percent + "%)"
u_patch = mpatches.Patch(color='darkgray', label=lab)
lab = "UIR - EE: " + str(uir_ee_num) + " (" + uir_ee_percent + "%), EN: " + str(uir_en_num) + " (" + uir_en_percent + "%), NE: " + str(uir_ne_num) + " (" + uir_ne_percent + "%), NN: " + str(uir_nn_num) + " (" + uir_nn_percent + "%)"
ur2_patch = mpatches.Patch(color='lightblue', label=lab)
plt.legend(handles=[di_patch, u_patch, ur2_patch], fontsize=10)
ax1.set_xlabel('Enrichment state pair tag')
ax1.set_ylabel('Read pair number')
plt.grid(True)
fig1.set_size_inches(7,5)
plt.savefig(pdf_name_boxplots_read_pair_numbers, bbox_inches='tight')
plt.close()


### Create barplots for proportions of interactions within EE, EN, NE and NN
############################################################################

def create_barplot_for_enrichment_pair_tag_percentages(title, i_percentages, i_numbers):

    fig, ax = plt.subplots()

    xticklables_labels = ['DI', 'UI', 'UIR']
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

ee_percentages = [round(100 * di_ee_num / di_num, 2), round(100 * ui_ee_num / ui_num, 2), round(100 * uir_ee_num / uir_num, 2)]
ee_numbers = [di_ee_num, ui_ee_num, uir_ee_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within EE", ee_percentages, ee_numbers))

en_percentages = [round(100 * di_en_num / di_num, 2), round(100 * ui_en_num / ui_num, 2), round(100 * uir_en_num / uir_num, 2)]
en_numbers = [di_en_num, ui_en_num, uir_en_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within EN", en_percentages, en_numbers))

ne_percentages = [round(100 * di_ne_num / di_num, 2), round(100 * ui_ne_num / ui_num, 2), round(100 * uir_ne_num / uir_num, 2)]
ne_numbers = [di_ne_num, ui_ne_num, uir_ne_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within NE", ne_percentages, ne_numbers))

nn_percentages = [round(100 * di_nn_num / di_num, 2), round(100 * ui_nn_num / ui_num, 2), round(100 * uir_nn_num / uir_num, 2)]
nn_numbers = [di_nn_num, ui_nn_num, uir_nn_num]
pdf.savefig(create_barplot_for_enrichment_pair_tag_percentages("Proportion of interactions within NN", nn_percentages, nn_numbers))

pdf.close()
