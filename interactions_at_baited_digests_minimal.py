"""
In this script the histograms with the paairwise differences of interaction distances of interactions at the same bait
are generated with minimal code.
"""

import argparse
import matplotlib.pyplot as plt
from diachr.diachromatic_interaction_set import DiachromaticInteractionSet

# Parse command line
####################

parser = argparse.ArgumentParser(description='Get BED files with regions spanned by interactions that end in baited digests.')
parser.add_argument('-o', '--out-prefix',
                    help='Common prefix for all generated files, which can also contain the path.',
                    default='OUT_PREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file',
                    help='Input file in Diachromatic11 interaction format.',
                    required=True)
parser.add_argument('-c', '--use-chromosome',
                    help='If specified, the analysis will be done for one chromosome only, e.g. chr15.',
                    default='USE_ALL')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.diachromatic_interaction_file
use_chromosome = args.use_chromosome

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --use-chromosome: " + use_chromosome)


class BaitedDigest:
    """
    Objects of this class are used to group interactions interactions that end in the same bait. Interactions are broken
    down by category and direction as seen from the bait. For the interaction categories directed (DI), undirected
    reference (UIR), undirected (UI) and all interactions (ALL) there are two lists of interactions each, one with
    interactions from the bait to the left (NE) and one with interactions from the bait to the right (EN).
    """

    # Dictionary with lists of interactions ending in baited digest
    interactions = None

    def __init__(self):

        self.interactions = {
            'DI': {
                'NE': [],
                'EN': []},
            'UIR': {
                'NE': [],
                'EN': []},
            'UI': {
                'NE': [],
                'EN': []},
            'ALL': {
                'NE': [],
                'EN': []}
        }

    def get_all_pairwise_differences_of_interaction_distances(self, i_cat, dir):

        d_inter_list = self.interactions[i_cat][dir]
        d_inter_list_len = len(d_inter_list)

        # Get list of all pairwise interaction distances
        pairwise_i_dist_diffs = []
        for i in range(0, d_inter_list_len):
            dist_a = d_inter_list[i].i_dist
            for j in range(i + 1, d_inter_list_len):
                dist_b = d_inter_list[j].i_dist
                pairwise_i_dist_diffs.append(abs(dist_a - dist_b))

        return pairwise_i_dist_diffs


# Read file with interactions
#############################

interaction_set = DiachromaticInteractionSet()
interaction_set.parse_file(diachromatic_interaction_file, verbose=True)


# Iterate interactions and assign them to BaitedDigest objects
##############################################################

print("[INFO] Iterating interactions and assigning them to BaitedDigest objects ...")

# Dictionary with BaitedDigest objects
baited_digests = dict()

# Variables for counting interactions
n_interaction = {'ALL_NN_EE': 0, 'ALL_NE_EN': 0, 'DI': 0, 'UIR':0, 'UI': 0}

# Iterate interactions
for d_inter in interaction_set.interaction_list:

    if use_chromosome != 'USE_ALL':
        if d_inter.chrA != use_chromosome:
            continue

    # We skip interactions that don't have an enriched digest
    if d_inter.enrichment_status_tag_pair == 'NN':
        n_interaction['ALL_NN_EE'] += 1
        continue

    # We also skip interactions that don't have two enriched digests
    if d_inter.enrichment_status_tag_pair == 'EE':
        n_interaction['ALL_NN_EE'] += 1
        continue

    # Count remaining NE and EN interactions
    n_interaction['ALL_NE_EN'] += 1
    n_interaction[d_inter.get_category()] += 1

    # Add interaction to list of BaitedDigest object
    if d_inter.enrichment_status_tag_pair == 'EN':
        coord_key_baited_digest = d_inter.chrA + '\t' + str(d_inter.fromA) + '\t' + str(d_inter.toA)

    if d_inter.enrichment_status_tag_pair == 'NE':
        coord_key_baited_digest = d_inter.chrB + '\t' + str(d_inter.fromB) + '\t' + str(d_inter.toB)

    if coord_key_baited_digest not in baited_digests:
        baited_digests[coord_key_baited_digest] = BaitedDigest()

    baited_digests[coord_key_baited_digest].interactions['ALL'][d_inter.enrichment_status_tag_pair].append(d_inter)
    baited_digests[coord_key_baited_digest].interactions[d_inter.get_category()][d_inter.enrichment_status_tag_pair].append(d_inter)

print("[INFO] ... done.")

print("[INFO] Interaction numbers")
print("\t[INFO] Total number of NN and EE interactions (not used): " + str(n_interaction['ALL_NN_EE']))
print("\t[INFO] Total number of NE and EN interactions (ALL): " + str(n_interaction['ALL_NE_EN']))
print("\t[INFO] Number of directed interactions (DI): " + str(n_interaction['DI']))
print("\t[INFO] Number of undirected reference interactions (UIR): " + str(n_interaction['UIR']))
print("\t[INFO] Number of undirected interactions (UI): " + str(n_interaction['UI']))


print("[INFO] Determining pairwise differences of interaction distances of interaction at the same bait ...")

# Create data structure
dist_diffs = dict()
i_cats = ['DI', 'UIR', 'UI', 'ALL']
enr_cats = ['NE', 'EN']
for i_cat in i_cats:
    dist_diffs[i_cat] = dict()
    for enr_cat in enr_cats:
        dist_diffs[i_cat][enr_cat] = []

# Get pairwise differences
for key in baited_digests:
    for i_cat in i_cats:
        for enr_cat in enr_cats:
            dist_diffs[i_cat][enr_cat].extend(
                baited_digests[key].get_all_pairwise_differences_of_interaction_distances(i_cat, enr_cat))

print("[INFO] ... done.")


# print("[INFO] Writing pairwise interaction distances to text files (for plotting in R) ...")
#
# for i_cat in i_cats:
#     for enr_cat in enr_cats:
#         output_stream = open(out_prefix + '_pairwise_dist_diff_' + i_cat + '_' + enr_cat + '.tab', 'wt')
#         for dist_diff in dist_diffs[i_cat][enr_cat]:
#             output_stream.write(str(dist_diff) + '\n')
#         output_stream.close()
#
# print("[INFO] ... done.")


print("[INFO] Creating histograms using matplotlib...")

# Prepare figure
fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(10, 12))
x_ticks=[0,500000,1000000,1500000,2000000]
bin_width = 20000
i_cat_colors = {'DI': (255/255,163/255,0/255,1), 'UIR': (171/255,215/255,230/255,1), 'UI': (210/255,210/255,210/255,1), 'ALL': 'black'}
i_cat_names = {'DI': 'Directed', 'UIR': 'Undirected reference', 'UI': 'Undirected', 'ALL': 'All'}

# Create histograms for DI, UIR, UI and ALL
for i in [0,1,2,3]:

    # Prepare bins
    x_max = max(dist_diffs[i_cats[i]]['NE'] + dist_diffs[i_cats[i]]['EN'])
    bins = range(0, x_max + bin_width, bin_width)

    # Create histogram for NE
    counts_1, bins, patches = ax[i][0].hist(
        dist_diffs[i_cats[i]]['NE'],
        bins=bins, density=False,
        facecolor=i_cat_colors[i_cats[i]],
        edgecolor="black", alpha=0.75)
    ax[i][0].set_title(i_cat_names[i_cats[i]] + ' to the left', loc='left')
    ax[i][0].set_xlabel('Pairwise differences of interaction distances')
    ax[i][0].set_ylabel('Frequency')
    ax[i][0].set_xticks(x_ticks)
    ax[i][0].set_xlim(0, 2000000)

    # Create histogram for NE
    counts_2, bins, patches = ax[i][1].hist(
        dist_diffs[i_cats[i]]['EN'],
        bins=bins, density=False,
        facecolor=i_cat_colors[i_cats[i]],
        edgecolor="black", alpha=0.75)
    ax[i][1].set_title(i_cat_names[i_cats[i]] + ' to the right', loc='left')
    ax[i][1].set_xlabel('Pairwise differences of interaction distances')
    ax[i][1].set_ylabel('Frequency')
    ax[i][1].set_xticks(x_ticks)
    ax[i][1].set_xlim(0, 2000000)

    # Make y-axes comparable
    ymax = max(max(counts_1),max(counts_2))
    ax[i][0].set_ylim(0, ymax)
    ax[i][1].set_ylim(0, ymax)

# Draw vertical lines for UIR - NE
ax[1][0].axvline(1 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(2 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(3 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(4 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(5 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(6 * 270500, color='blue', linewidth=0.5, zorder=0)
ax[1][0].axvline(7 * 270500, color='blue', linewidth=0.5, zorder=0)

# Save figure
fig.tight_layout()
fig.savefig(out_prefix + "_pairwise_difference_historgrams.pdf")

print("[INFO] ... done.")
