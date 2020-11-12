import os
import argparse
from diachr import BinomialInteractionModel
from diachr import EnhancedInteraction
from diachr import EnhancedInteractionParser
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Explore binomial model using simuated and real interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--i-num', help='Number of simulated interactions', default=10000)
parser.add_argument('--n-max', help='Simulate interactions with 1 to n read pairs.',default=500)
parser.add_argument('--p-value-cutoff', help='P-value threshold used for the classification of directed and undirected interactions.', default=0.05)
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
n_max = int(args.n_max)
i_num = int(args.i_num)
p_value_cutoff = float(args.p_value_cutoff)
enhanced_interaction_file = args.enhanced_interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --i-num: " + str(i_num))
print("\t[INFO] --n-max: " + str(n_max))
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))
print("\t[INFO] --enhanced-interaction-file: " + str(enhanced_interaction_file))

# Create object
bim = BinomialInteractionModel(p_value_cutoff)

# Simulate interactions
n_max=400
i_num=100000000
bim._p_value_cutoff = 0.05
n_list_s, n_sim_list, n_sig_list = bim.simulate_interactions(n_max=n_max, i_num=i_num)

plt.rcParams["figure.figsize"] = (10,5)
fig, ax1 = plt.subplots()
color = 'gray'
ax1.set_xlabel('Number of read pairs per interaction (n)')
ax1.set_ylabel('Simulated interactions', color=color)
ax1.scatter(n_list_s, n_sim_list, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = (255/255,163/255,0/255,1)
ax2.set_ylabel('Simulated significant interactions', color=color)  # we already handled the x-label with ax1
ax2.plot(n_list_s, n_sig_list, color=color, linewidth=0.5)
ax2.scatter(n_list_s, n_sig_list, color=color, alpha=0.5)
ax2.hlines(bim._p_value_cutoff*(i_num/n_max),0, n_max, color='red', linestyle='dashed')
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(bim._out_prefix + "_simulated_interactions.pdf", format = "pdf")

# plt.scatter(n_list_s, n_sim_list, alpha=0.5)
# plt.title('Simulated interactions')
# plt.xlabel('Number of read pairs per interaction (n)')
# plt.ylabel('Number of interactions')
#
# # Scatterplot
# plt.scatter(n_list_s, n_sig_list, alpha=0.5)
# plt.title('Simulated significant interactions')
# plt.xlabel('Number of read pairs per interaction (n)')
# plt.ylabel('Number of interactions')
#
# plt.show()

# print("[INFO] " + "Generating plot: Significant simulated interactions versus n ...")
# plt.plot(signum_list)
# plt.grid(True)
# plt.xlabel("n")
# plt.ylabel("# Significant interactions")
# sub_title = "# Interactions: " + str(i_num) + " | Max n: " + str(n_max) + " | P-value cutoff: " + str(p_value_cutoff)
# plt.suptitle(sub_title)
# plt.savefig(pdf_name, format = "pdf")


# Count DI, UIR and UI for each n in enhanced interaction file
if enhanced_interaction_file != None:
    if not os.path.exists(enhanced_interaction_file):
        raise FileNotFoundError("Could not find IE file")
    n_list_e, n_di_list, n_uir_list, n_ui_list = bim.count_di_uir_and_ui_for_each_n(ei_file=enhanced_interaction_file)




# bim.write_simulated_interaction_counts(outprefix=out_prefix)
#
# bim.write_simulated_interaction_counts()
# if diachromatic_interaction_file != None:
#     print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)
#     bim.write_significant_empirical_interactions(ei_file=diachromatic_interaction_file, n_indef=n_max)

print("******** Done new version **********")