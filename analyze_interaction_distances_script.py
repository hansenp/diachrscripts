#!/usr/bin/env python
import diachrscripts_toolkit as dclass
import argparse
import gzip
import statistics as statistics
import matplotlib.pyplot as plt
import numpy as np

# Parse command line arguments
parser = argparse.ArgumentParser(description='Analyze distances between interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--status-pair-flag', help='Pair of \'A\' and \'I\' depending on whether a digest was selected for enrichment (A) or not (I).')
parser.add_argument('--p-value-cutoff', help='P-value cutoff used for categorization of interactions.')


args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
status_pair_flag = args.status_pair_flag
if status_pair_flag != "ALL":
    status_pair_flag = sorted(status_pair_flag)[0] + sorted(status_pair_flag)[1]
p_value_cutoff = float(args.p_value_cutoff)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --interaction-file: " + diachromatic_interaction_file)
print("\t[INFO] --status-pair-flag: " + status_pair_flag)
print("\t[INFO] --p-value-cutoff: " + str(p_value_cutoff))

# Init arrays for distribution of distances between digests
distance_array_simple = []
distance_array_twisted = []
distance_array_undirected = []
distance_array_indefinable = []

n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_indefinable_interaction = 0

n_simple_interaction_rp = 0
n_twisted_interaction_rp = 0
n_undirected_interaction_rp = 0
n_indefinable_interaction_rp = 0

n_trans_short_range_interaction = 0
n_non_status_pair_flag_interaction = 0

# Iterate interactions and collect distances
print("[INFO] Determining distances between digest pairs in " + diachromatic_interaction_file + " ...")
n_interaction_total = 0
with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
    line = fp.readline()
    while line:

        if n_interaction_total%1000000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")
        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)

        # restrict analysis to interactions between targeted promoters
        if status_pair_flag != "ALL" and interaction.get_digest_status_pair_flag() != status_pair_flag:
            n_non_status_pair_flag_interaction +=1
            line = fp.readline()
            continue

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(10000)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        # Determine distance between the centers of digests
        distance = interaction.get_digest_distance()

        # set the type of interaction based on P-value ('S', 'T', 'U', 'NA')
        if interaction.get_interaction_type() == "TBD":
            interaction.set_interaction_type("TBD", p_value_cutoff)

        n_read_pair = interaction.n_simple + interaction.n_twisted

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "NA":
            distance_array_indefinable.append(distance)
            n_indefinable_interaction += 1
            n_indefinable_interaction_rp = n_indefinable_interaction_rp + n_read_pair
        elif interaction.get_interaction_type() == "U":
            distance_array_undirected.append(distance)
            n_undirected_interaction += 1
            n_undirected_interaction_rp = n_undirected_interaction_rp + n_read_pair
        elif interaction.get_interaction_type() == "S":
            distance_array_simple.append(distance)
            n_simple_interaction += 1
            n_simple_interaction_rp = n_simple_interaction_rp + n_read_pair
        elif interaction.get_interaction_type() == "T":
            distance_array_twisted.append(distance)
            n_twisted_interaction += 1
            n_twisted_interaction_rp = n_twisted_interaction_rp + n_read_pair
        else:
            line = fp.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        line = fp.readline()

fp.close()

# Output summary
print("[INFO] " + "Summary statistics")
print("\t[INFO] Total number of interactions: " + str(n_interaction_total))
print("\t[INFO] Number of non " + status_pair_flag + " interactions: " + str(n_non_status_pair_flag_interaction) + " (discarded)")
print("\t[INFO] Number of trans and short range interactions: " + str(n_trans_short_range_interaction) + " (discarded)")
print("\t[INFO] Number of simple interactions: " + str(n_simple_interaction))
print("\t[INFO] Number of twisted interactions: " + str(n_twisted_interaction))
print("\t[INFO] Number of undirected interactions: " + str(n_undirected_interaction))
print("\t[INFO] Number of indefinable interactions: " + str(n_indefinable_interaction))
print("\t[INFO] Number of read pairs in simple interactions: " + str(n_simple_interaction_rp) + " (" + str(n_simple_interaction_rp/n_simple_interaction) + ")")
print("\t[INFO] Number of read pairs in twisted interactions: " + str(n_twisted_interaction_rp) + " (" + str(n_twisted_interaction_rp/n_twisted_interaction) + ")")
print("\t[INFO] Number of read pairs in undirected interactions: " + str(n_undirected_interaction_rp) + " (" + str(n_undirected_interaction_rp/n_undirected_interaction) + ")")
print("\t[INFO] Number of read pairs in indefinable interactions: " + str(n_indefinable_interaction_rp) + " (" + str(n_indefinable_interaction_rp/n_indefinable_interaction) + ")")

print("[INFO] " + "Writing numpy arrays with distances to disk ...")
file_path_name = out_prefix + "_distance_array_simple"
np.save(file_path_name, np.array(distance_array_simple))
file_path_name = out_prefix + "_distance_array_twisted"
np.save(file_path_name, np.array(distance_array_twisted))
file_path_name = out_prefix + "_distance_array_undirected"
np.save(file_path_name, np.array(distance_array_undirected))
file_path_name = out_prefix + "_distance_array_indefinable"
np.save(file_path_name , np.array(distance_array_indefinable))

print("[INFO] " + "Writing arrays with distances to files ...")
file_name = out_prefix + "_digest_distances_simple.txt"
f_output = open(file_name, 'wt')
for d in distance_array_simple:
    f_output.write(str(d) + "\n")
f_output.close()
file_name = out_prefix + "_digest_distances_twisted.txt"
f_output = open(file_name, 'wt')
for d in distance_array_twisted:
    f_output.write(str(d) + "\n")
f_output.close()
file_name = out_prefix + "_digest_distances_undirected.txt"
f_output = open(file_name, 'wt')
for d in distance_array_undirected:
    f_output.write(str(d) + "\n")
f_output.close()
file_name = out_prefix + "_digest_distances_indefinable.txt"
f_output = open(file_name, 'wt')
for d in distance_array_indefinable:
    f_output.write(str(d) + "\n")
f_output.close()

print("[INFO] " + "... done.")
