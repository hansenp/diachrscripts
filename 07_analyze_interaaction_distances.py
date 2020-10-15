#!/usr/bin/env python

"""
This script takes an enhanced interaction (EI) file (created with the script '06_xxxxxxxxxxxxxxxxxxxxxxxxxx') with the
following interaction categories:

   1. DI - Directed interactions
   2. UI - Undirected interactions

that are specified in the third column.

Based on the enrichment pair tag in column X the distribution of interaction distances will be determined for

Directed interactions: EE, EN, NE and NN
Unndirected interactions: EE, EN, NE and NN


Finally, summary statistics about numbers of interactions within different enrichment state categories (EE, AI and II)
and corresponding read pair numbers are printed to the screen and two boxplots are created, one for read pair numbers
and another one for interaction distances:

   '<OUT_PREFIX>_read_pair_number_boxplot.pdf'
   '<OUT_PREFIX>_interaction_distance_boxplot.pdf'

Furthermore, barplots for the proportions of interactions and associated digests are created:

   '<OUT_PREFIX>_interaction_enrichment_pair_tags_barplot.pdf'
   '<OUT_PREFIX>__digest_enrichment_tags_barplot.pdf'

"""


### Create boxplots for interaction distances
#############################################

plt.rcParams.update({'font.size':15})

data = [
    dir_inter_aa_dist_array, dir_inter_ai_dist_array, dir_inter_ii_dist_array,
    undir_inter_aa_dist_array,undir_inter_ai_dist_array, undir_inter_ii_dist_array,
    undir_ref_2_inter_aa_dist_array, undir_ref_2_inter_ai_dist_array, undir_ref_2_inter_ii_dist_array
]
labels = ['EE', 'EN', 'NN', 'EE', 'EN', 'NN','EE', 'EN', 'NN']
fig1, ax1 = plt.subplots()
ax1.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax1.set_title('Distributions of interaction distances for DI, U, and UR (' + out_prefix + ')')
ax1.boxplot(data, showfliers=False, labels=labels)

box = ax1.boxplot(data, showfliers=False, labels=labels, patch_artist=True)

colors = ['orange', 'orange', 'orange', 'darkgray', 'darkgray', 'darkgray', 'lightblue', 'lightblue', 'lightblue']

for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

lab = "DI - AA: " + str(dir_inter_aa_num) + " (" + dir_inter_aa_percentage + "%), AI: " + str(dir_inter_ai_num) + " (" + dir_inter_ai_percentage + "%), II: " + str(dir_inter_ii_num) + " (" + dir_inter_ii_percentage + "%)"
di_patch = mpatches.Patch(color='orange', label=lab)
lab = "U - AA: " + str(undir_inter_aa_num) + " (" + undir_inter_aa_percentage + "%), AI: " + str(undir_inter_ai_num) + " (" + undir_inter_ai_percentage + "%), II: " + str(undir_inter_ii_num) + " (" + undir_inter_ii_percentage + "%)"
u_patch = mpatches.Patch(color='darkgray', label=lab)
lab = "UR 2 - AA: " + str(undir_ref_2_inter_aa_num) + " (" + undir_ref_2_inter_aa_percentage + "%), AI: " + str(undir_ref_2_inter_ai_num) + " (" + undir_ref_2_inter_ai_percentage + "%), II: " + str(undir_ref_2_inter_ii_num) + " (" + undir_ref_2_inter_ii_percentage + "%)"
ur2_patch = mpatches.Patch(color='lightblue', label=lab)
plt.legend(handles=[di_patch, u_patch, ur2_patch])
ax1.set_xlabel('Enrichment state pair tag')
ax1.set_ylabel('Interaction distance')
plt.grid(True)
fig1.set_size_inches(7,5)
plt.savefig(pdf_name_boxplots_interaction_distances, bbox_inches='tight')
plt.close()
