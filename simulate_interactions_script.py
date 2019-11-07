import argparse
import gzip
import diachrscripts_toolkit as dclass
import numpy as np

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine expression category levels for interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file

print("[INFO] " + "Input parameters")
print("\t[INFO] Analysis for: " + out_prefix)
print("\t[INFO] Interaction file: " + diachromatic_interaction_file)

### Start execution
###################

n_interaction_total = 0
n_trans_short_range_interaction = 0

n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_indefinable_interaction = 0

n_simple_interaction_sim = 0
n_twisted_interaction_sim = 0
n_undirected_interaction_sim = 0
n_indefinable_interaction_sim = 0

# Init arrays for qq plot of P-values
pvals_original = []
pvals_simulated = []

# open two streams for gzipped output of original and simulated reads
file_name_original = out_prefix + "_original_interactions.tsv.gz"
f_output_original = gzip.open(file_name_original, 'wt')
file_name_simulated = out_prefix + "_simulated_interactions.tsv.gz"
f_output_simulated = gzip.open(file_name_simulated, 'wt')

# iterate interactions
print("[INFO] Simulating random simple twisted read pairs for each interaction in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, 'r' + 't') as f_input:

    line = f_input.readline()

    while line:

        if n_interaction_total%100000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)
        sim_interaction = interaction.get_simulated_copy()

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(1000)):
            n_trans_short_range_interaction += 1
            line = f_input.readline()
            continue

        # set the type of interaction based on P-value ('S', 'T', 'U', 'NA')
        if interaction.get_interaction_type() == "TBD":
            interaction.set_interaction_type("TBD", 0.05)

        if sim_interaction.get_interaction_type() == "TBD":
            sim_interaction.set_interaction_type("TBD", 0.05)

        # count interaction type for original interaction
        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "NA":
            n_indefinable_interaction += 1
        elif interaction.get_interaction_type() == "U":
            n_undirected_interaction += 1
        elif interaction.get_interaction_type() == "S":
            n_simple_interaction += 1
        elif interaction.get_interaction_type() == "T":
            n_twisted_interaction += 1
        else:
            line = f_input.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        # count interaction type for simulated interaction
        if sim_interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif sim_interaction.get_interaction_type() == "NA":
            n_indefinable_interaction_sim += 1
        elif sim_interaction.get_interaction_type() == "U":
            n_undirected_interaction_sim += 1
        elif sim_interaction.get_interaction_type() == "S":
            n_simple_interaction_sim += 1
        elif sim_interaction.get_interaction_type() == "T":
            n_twisted_interaction_sim += 1
        else:
            line = f_input.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        # write interactions to file
        interaction_coords = interaction.get_first_digest().chromosome + '\t'\
                                + str(interaction.get_first_digest().start) + '\t'\
                                + str(interaction.get_first_digest().end) + '\t' \
                                + interaction.get_second_digest().chromosome + '\t' \
                                + str(interaction.get_second_digest().start) + '\t' \
                                + str(interaction.get_second_digest().end) + '\t'
        f_output_original.write(interaction_coords + '\t' + str(interaction.n_simple) + '\t' + str(interaction.n_twisted) + '\t' + str(interaction.get_binomial_p_value()) + '\n')
        f_output_simulated.write(interaction_coords + '\t' + str(sim_interaction.n_simple) + '\t' + str(sim_interaction.n_twisted) + '\t' + str(sim_interaction.get_binomial_p_value()) + '\n')

        # add P-values to array for qq plot
        pvals_original.append(interaction.get_binomial_p_value())
        pvals_simulated.append(sim_interaction.get_binomial_p_value())

        line = f_input.readline()

    print("... done.")

f_input.close()
f_output_original.close()
f_output_simulated.close()

# output summary
print("[INFO] " + "Summary statistics")
print("\t[INFO] Total number of interactions: " + str(n_interaction_total))
print("\t[INFO] Number of trans short range interactions (discarded): " + str(n_trans_short_range_interaction))
print("\t[INFO] Number of simple interactions (original vs. simulated): " + str(n_simple_interaction) + "\t" + str(n_simple_interaction_sim))
print("\t[INFO] Number of twisted interactions (original vs. simulated): " + str(n_twisted_interaction) + "\t" + str(n_twisted_interaction_sim))
print("\t[INFO] Number of undirected interactions (original vs. simulated): " + str(n_undirected_interaction) + "\t" + str(n_undirected_interaction_sim))
print("\t[INFO] Number of indefinable interactions (original vs. simulated): " + str(n_indefinable_interaction) + "\t" + str(n_indefinable_interaction_sim))

# save arrays with P-values to disk so as we can use them in a notebook
print("[INFO] " + "Writing numpy arrays with distances to disk ...")
file_path_name_original = out_prefix + "_original_pvals"
np.save(file_path_name_original, np.array(pvals_original))
file_path_name_simulated = out_prefix + "_simulated_pvals"
np.save(file_path_name_simulated, np.array(pvals_simulated))
print("[INFO] " + "... done.")

print("[INFO] " + "Two files with interactions were generated:")
print("\t[INFO] Original interactions with P-values: " + file_name_original)
print("\t[INFO] Simulated interactions with P-values: " + file_name_simulated)

print("[INFO] " + "Two numpy arrays with P-values were saved to disk:")
print("\t[INFO] P-values of original data: " + file_path_name_original)
print("\t[INFO] P-values of simulated data: " + file_path_name_simulated)
