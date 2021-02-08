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

Since, in the EI format, the two digests of each interaction are in sequential order, interaactions with the tag NE go
to the left from the bait aand interaactions with the tag EN go to the right.

This script iterates the input EI file and collects the interaction distances for DI, UI and UIR separately for EE, NE,
EN and NN. The twelve resulting distributions are then shown in a boxplot with the following name:

   '<OUT_PREFIX>_interaction_distance_boxplot.pdf'

"""
####X The section above explains the input for the script, what it does and what files will be created.


import argparse
import gzip
import diachrscripts_toolkit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

### Parse command line
######################

parser = argparse.ArgumentParser(description='Select reference interactions from exclusive undirected interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file created with \'05_define_di_uie_and_uii.py\'. Column 3 contains either \'DI\' or\'UIE\'.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + enhanced_interaction_file)


### Prepare variables
#####################
####X The essential variables are declared in the following section.

# Directed interactions
# ---------------------

# Total number of interactions
di_num = 0

# Numbers of interactions within EE, NE, EN and NN
di_ee_num = 0
di_ne_num = 0
di_en_num = 0
di_nn_num = 0

# Arrays for interaction distances within EE, NE, EN and NN
di_ee_dist_array = []
di_ne_dist_array = []
di_en_dist_array = []
di_nn_dist_array = []

# Arrays for interaction distances at the level of simple and twisted read pairs within EE, NE, EN and NN
## Simple
di_ee_s_rp_dist_array = []
di_ne_s_rp_dist_array = []
di_en_s_rp_dist_array = []
di_nn_s_rp_dist_array = []
## Twisted
di_ee_t_rp_dist_array = []
di_ne_t_rp_dist_array = []
di_en_t_rp_dist_array = []
di_nn_t_rp_dist_array = []

# Directed simple interactions (is a subset of all directed interactions)
# -----------------------------------------------------------------------

# Total number of interactions
di_s_num = 0

# Numbers of interactions within EE, NE, EN and NN
di_ee_s_num = 0
di_ne_s_num = 0
di_en_s_num = 0
di_nn_s_num = 0

# Arrays for interaction distances within EE, NE, EN and NN
di_ee_s_dist_array = []
di_ne_s_dist_array = []
di_en_s_dist_array = []
di_nn_s_dist_array = []

# Directed simple interactions (is a subset of all directed interactions)
# -----------------------------------------------------------------------

# Total number of interactions
di_t_num = 0

# Numbers of interactions within EE, NE, EN and NN
di_ee_t_num = 0
di_ne_t_num = 0
di_en_t_num = 0
di_nn_t_num = 0

# Arrays for interaction distances within EE, NE, EN and NN
di_ee_t_dist_array = []
di_ne_t_dist_array = []
di_en_t_dist_array = []
di_nn_t_dist_array = []

# Undirected reference interactions
# ---------------------------------

# Total number of interactions
uir_num = 0

# Numbers of interactions within EE, NE, EN and NN
uir_ee_num = 0
uir_ne_num = 0
uir_en_num = 0
uir_nn_num = 0

# Arrays for interaction distances within EE, NE, EN and NN
uir_ee_dist_array = []
uir_ne_dist_array = []
uir_en_dist_array = []
uir_nn_dist_array = []

# Arrays for interaction distances at the level of simple and twisted read pairs within EE, NE, EN and NN
## Simple
uir_ee_s_rp_dist_array = []
uir_ne_s_rp_dist_array = []
uir_en_s_rp_dist_array = []
uir_nn_s_rp_dist_array = []
## Twisted
uir_ee_t_rp_dist_array = []
uir_ne_t_rp_dist_array = []
uir_en_t_rp_dist_array = []
uir_nn_t_rp_dist_array = []

# Undirected interactions
# -----------------------

# Total number of interactions
ui_num = 0

# Numbers of interactions within EE, NE, EN and NN
ui_ee_num = 0
ui_ne_num = 0
ui_en_num = 0
ui_nn_num = 0

# Arrays for interaction distances within EE, NE, EN and NN
ui_ee_dist_array = []
ui_ne_dist_array = []
ui_en_dist_array = []
ui_nn_dist_array = []

# Arrays for interaction distances at the level of simple and twisted read pairs within EE, NE, EN and NN
## Simple
ui_ee_s_rp_dist_array = []
ui_ne_s_rp_dist_array = []
ui_en_s_rp_dist_array = []
ui_nn_s_rp_dist_array = []
## Twisted
ui_ee_t_rp_dist_array = []
ui_ne_t_rp_dist_array = []
ui_en_t_rp_dist_array = []
ui_nn_t_rp_dist_array = []

### Prepare output files
########################

# PDF file with boxplots for distributions of interaction distances
pdf_name_boxplots_read_pair_numbers = out_prefix + "_interaction_distances_boxplot.pdf"


### 1st pass: Determine distribution of interaction distances for DI
####################################################################

print("[INFO] Collecting information about distances ...")

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

        n_simple = int(line.split('\t')[4].split(':')[0])
        n_twisted = int(line.split('\t')[4].split(':')[1])

        # Collect read pair numbers for DI, UI and UIR
        if interaction_category == 'DI':

            di_num += 1
            if n_twisted < n_simple:
                di_s_num += 1
            else:
                di_t_num += 1

            if enrichment_pair_tag == 'EE':
                di_ee_num +=1
                di_ee_dist_array.append(i_dist)

                if n_twisted < n_simple:
                    di_ee_s_num += 1
                    di_ee_s_dist_array.append(i_dist)
                else:
                    di_ee_t_num += 1
                    di_ee_t_dist_array.append(i_dist)

                for s in range(0,n_simple):
                    di_ee_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    di_ee_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NE':
                di_ne_num += 1
                di_ne_dist_array.append(i_dist)

                if n_twisted < n_simple:
                    di_ne_s_num += 1
                    di_ne_s_dist_array.append(i_dist)
                else:
                    di_ne_t_num += 1
                    di_ne_t_dist_array.append(i_dist)

                for s in range(0,n_simple):
                    di_ne_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    di_ne_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'EN':
                di_en_num += 1
                di_en_dist_array.append(i_dist)

                if n_twisted < n_simple:
                    di_en_s_num += 1
                    di_en_s_dist_array.append(i_dist)
                else:
                    di_en_t_num += 1
                    di_en_t_dist_array.append(i_dist)

                for s in range(0,n_simple):
                    di_en_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    di_en_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NN':
                di_nn_num += 1
                di_nn_dist_array.append(i_dist)

                if n_twisted < n_simple:
                    di_nn_s_num += 1
                    di_nn_s_dist_array.append(i_dist)
                else:
                    di_nn_t_num += 1
                    di_nn_t_dist_array.append(i_dist)

                for s in range(0,n_simple):
                    di_nn_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    di_nn_t_rp_dist_array.append(i_dist)

        elif interaction_category == 'UIR':

            uir_num += 1

            if enrichment_pair_tag == 'EE':
                uir_ee_num += 1
                uir_ee_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    uir_ee_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    uir_ee_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NE':
                uir_ne_num += 1
                uir_ne_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    uir_ne_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    uir_ne_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'EN':
                uir_en_num += 1
                uir_en_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    uir_en_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    uir_en_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NN':
                uir_nn_num += 1
                uir_nn_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    uir_nn_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    uir_nn_t_rp_dist_array.append(i_dist)

        elif interaction_category == 'UI':

            ui_num += 1

            if enrichment_pair_tag == 'EE':
                ui_ee_num += 1
                ui_ee_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    ui_ee_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    ui_ee_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NE':
                ui_ne_num += 1
                ui_ne_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    ui_ne_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    ui_ne_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'EN':
                ui_en_num += 1
                ui_en_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    ui_en_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    ui_en_t_rp_dist_array.append(i_dist)

            elif enrichment_pair_tag == 'NN':
                ui_nn_num += 1
                ui_nn_dist_array.append(i_dist)
                for s in range(0,n_simple):
                    ui_nn_s_rp_dist_array.append(i_dist)
                for t in range(0,n_twisted):
                    ui_nn_t_rp_dist_array.append(i_dist)

        else:
            print("[Error] Interaction category must be either \'DI\', \'UI\' or \'UIR\'!")
            exit(1)

        line = fp.readline()

print("\t[INFO] ... done.")

### Create boxplots for interaction distances
#############################################

di_ee_percent = "{0:.2f}".format(100 * di_ee_num / di_num)
di_en_percent = "{0:.2f}".format(100 * di_en_num / di_num)
di_ne_percent = "{0:.2f}".format(100 * di_ne_num / di_num)
di_nn_percent = "{0:.2f}".format(100 * di_nn_num / di_num)

ui_ee_percent = "{0:.2f}".format(100 * ui_ee_num / ui_num)
ui_en_percent = "{0:.2f}".format(100 * ui_en_num / ui_num)
ui_ne_percent = "{0:.2f}".format(100 * ui_ne_num / ui_num)
ui_nn_percent = "{0:.2f}".format(100 * ui_nn_num / ui_num)

uir_ee_percent = "{0:.2f}".format(100 * uir_ee_num / uir_num)
uir_en_percent = "{0:.2f}".format(100 * uir_en_num / uir_num)
uir_ne_percent = "{0:.2f}".format(100 * uir_ne_num / uir_num)
uir_nn_percent = "{0:.2f}".format(100 * uir_nn_num / uir_num)

plt.rcParams.update({'font.size':15})

data = [
    di_ee_dist_array, di_en_dist_array, di_ne_dist_array, di_nn_dist_array,
    ui_ee_dist_array, ui_en_dist_array, ui_ne_dist_array, ui_nn_dist_array,
    uir_ee_dist_array, uir_en_dist_array, uir_ne_dist_array, uir_nn_dist_array
]
labels = ['EE', 'EN', 'NE', 'NN', 'EE', 'EN', 'NE', 'NN','EE', 'EN', 'NE', 'NN']
fig1, ax1 = plt.subplots()
ax1.set_title('Distributions of interaction distances for DI, UI, and UIR (' + out_prefix + ')', fontsize=10)
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
ax1.set_ylabel('Interaction distance')
plt.grid(True)
fig1.set_size_inches(7,5)
plt.savefig(pdf_name_boxplots_read_pair_numbers, bbox_inches='tight')
plt.close()


# Write interaction distances to text files for downstream analyzes
# -----------------------------------------------------------------

print("[INFO] Writing interaction distances to files ...")

def wrtite_interaction_distances_to_file(distance_array, file_name):
    tab_stream_output = open(file_name, 'wt')
    for d in distance_array:
        tab_stream_output.write(str(d) + '\n')
    tab_stream_output.close()


print("\t[INFO] Directed interactions")
wrtite_interaction_distances_to_file(di_ee_dist_array, out_prefix + "_di_ee_dist_array.tab")
wrtite_interaction_distances_to_file(di_ne_dist_array, out_prefix + "_di_ne_dist_array.tab")
wrtite_interaction_distances_to_file(di_en_dist_array, out_prefix + "_di_en_dist_array.tab")
wrtite_interaction_distances_to_file(di_nn_dist_array, out_prefix + "_di_nn_dist_array.tab")

print("\t[INFO] Directed simple interactions")
wrtite_interaction_distances_to_file(di_ee_s_dist_array, out_prefix + "_di_ee_s_dist_array.tab")
wrtite_interaction_distances_to_file(di_ne_s_dist_array, out_prefix + "_di_ne_s_dist_array.tab")
wrtite_interaction_distances_to_file(di_en_s_dist_array, out_prefix + "_di_en_s_dist_array.tab")
wrtite_interaction_distances_to_file(di_nn_s_dist_array, out_prefix + "_di_nn_s_dist_array.tab")

print("\t[INFO] Directed twisted interactions")
wrtite_interaction_distances_to_file(di_ee_t_dist_array, out_prefix + "_di_ee_t_dist_array.tab")
wrtite_interaction_distances_to_file(di_ne_t_dist_array, out_prefix + "_di_ne_t_dist_array.tab")
wrtite_interaction_distances_to_file(di_en_t_dist_array, out_prefix + "_di_en_t_dist_array.tab")
wrtite_interaction_distances_to_file(di_nn_t_dist_array, out_prefix + "_di_nn_t_dist_array.tab")

print("\t[INFO] Undirected reference interactions")
wrtite_interaction_distances_to_file(uir_ee_dist_array, out_prefix + "_uir_ee_dist_array.tab")
wrtite_interaction_distances_to_file(uir_ne_dist_array, out_prefix + "_uir_ne_dist_array.tab")
wrtite_interaction_distances_to_file(uir_en_dist_array, out_prefix + "_uir_en_dist_array.tab")
wrtite_interaction_distances_to_file(uir_nn_dist_array, out_prefix + "_uir_nn_dist_array.tab")

print("\t[INFO] Undirected interactions")
wrtite_interaction_distances_to_file(ui_ee_dist_array, out_prefix + "_ui_ee_dist_array.tab")
wrtite_interaction_distances_to_file(ui_ne_dist_array, out_prefix + "_ui_ne_dist_array.tab")
wrtite_interaction_distances_to_file(ui_en_dist_array, out_prefix + "_ui_en_dist_array.tab")
wrtite_interaction_distances_to_file(ui_nn_dist_array, out_prefix + "_ui_nn_dist_array.tab")

###############################################################################
# Interaction distances at the level of read pairs

print("\t[INFO] Directed interactions - Simple read pairs")
wrtite_interaction_distances_to_file(di_ee_s_rp_dist_array, out_prefix + "_di_ee_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_ne_s_rp_dist_array, out_prefix + "_di_ne_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_en_s_rp_dist_array, out_prefix + "_di_en_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_nn_s_rp_dist_array, out_prefix + "_di_nn_s_rp_dist_array.tab")

print("\t[INFO] Directed interactions - Twisted read pairs")
wrtite_interaction_distances_to_file(di_ee_t_rp_dist_array, out_prefix + "_di_ee_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_ne_t_rp_dist_array, out_prefix + "_di_ne_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_en_t_rp_dist_array, out_prefix + "_di_en_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(di_nn_t_rp_dist_array, out_prefix + "_di_nn_t_rp_dist_array.tab")

print("\t[INFO] Undirected reference interactions - Simple read pairs")
wrtite_interaction_distances_to_file(uir_ee_s_rp_dist_array, out_prefix + "_uir_ee_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_ne_s_rp_dist_array, out_prefix + "_uir_ne_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_en_s_rp_dist_array, out_prefix + "_uir_en_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_nn_s_rp_dist_array, out_prefix + "_uir_nn_s_rp_dist_array.tab")

print("\t[INFO] Undirected reference interactions - Twisted read pairs")
wrtite_interaction_distances_to_file(uir_ee_t_rp_dist_array, out_prefix + "_uir_ee_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_ne_t_rp_dist_array, out_prefix + "_uir_ne_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_en_t_rp_dist_array, out_prefix + "_uir_en_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(uir_nn_t_rp_dist_array, out_prefix + "_uir_nn_t_rp_dist_array.tab")

print("\t[INFO] Undirected interactions - Simple read pairs")
wrtite_interaction_distances_to_file(ui_ee_s_rp_dist_array, out_prefix + "_ui_ee_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_ne_s_rp_dist_array, out_prefix + "_ui_ne_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_en_s_rp_dist_array, out_prefix + "_ui_en_s_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_nn_s_rp_dist_array, out_prefix + "_ui_nn_s_rp_dist_array.tab")

print("\t[INFO] Undirected interactions - Twisted read pairs")
wrtite_interaction_distances_to_file(ui_ee_t_rp_dist_array, out_prefix + "_ui_ee_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_ne_t_rp_dist_array, out_prefix + "_ui_ne_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_en_t_rp_dist_array, out_prefix + "_ui_en_t_rp_dist_array.tab")
wrtite_interaction_distances_to_file(ui_nn_t_rp_dist_array, out_prefix + "_ui_nn_t_rp_dist_array.tab")

print("[INFO] ... done.")
