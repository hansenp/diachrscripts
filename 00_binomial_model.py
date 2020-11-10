import os
import argparse
from diachr import BinomialInteractionModel
from diachr import EnhancedInteraction
from diachr import EnhancedInteractionParser

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


# bim = BinomialInteractionModel(n_max=n_max, 
#             i_num=i_num, 
#             p_value_cutoff=p_value_cutoff, 
#             out_prefix=out_prefix)

bim_2 = BinomialInteractionModel(p_value_cutoff)
n_list, n_sim_list, n_sig_list = bim_2.count_simulated_interactions(n_max=n_max, i_num=i_num)

# This should be an enhanced interaction file
if not os.path.exists(enhanced_interaction_file):
    raise FileNotFoundError("Could not find IE file")

n_list, n_di_list, n_uir_list, n_ui_list = bim_2.count_significant_empirical_interactions(ei_file=enhanced_interaction_file)



# bim.write_simulated_interaction_counts(outprefix=out_prefix)
#
# bim.write_simulated_interaction_counts()
# if diachromatic_interaction_file != None:
#     print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)
#     bim.write_significant_empirical_interactions(ei_file=diachromatic_interaction_file, n_indef=n_max)

print("******** Done new version **********")