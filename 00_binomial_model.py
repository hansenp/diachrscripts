import argparse
from diachr import BinomialInteractionModel

parser = argparse.ArgumentParser(description='Explore binomial model using simuated and real interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--i-num', help='Number of simulated interactions', default=10000)
parser.add_argument('--n-max', help='Simulate interactions with 1 to n read pairs.',default=500)
parser.add_argument('--p-value-cutoff', help='P-value threshold used for the classification of directed and undirected interactions.', default=0.05)
parser.add_argument('-d','--diachromatic-interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
n_max = int(args.n_max)
i_num = int(args.i_num)
p_value_cutoff = float(args.p_value_cutoff)
diachromatic_interaction_file = args.diachromatic_interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --i-num: " + str(i_num))
print("\t[INFO] --n-max: " + str(n_max))
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))
if diachromatic_interaction_file != None:
    print("\t[INFO] --diachromatic-interaction-file: " + diachromatic_interaction_file)

bim = BinomialInteractionModel(n_max=n_max, 
            i_num=i_num, 
            p_value_cutoff=p_value_cutoff, 
            out_prefix=out_prefix, 
            diachromatic_interaction_file=diachromatic_interaction_file)

bim.do_simulation()

print("******** Done new version **********")