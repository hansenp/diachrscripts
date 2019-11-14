#!/usr/bin/env python
import diachrscripts_toolkit as dclass
import argparse
import gzip
import statistics as statistics
import matplotlib.pyplot as plt
import numpy as np

### Parse command line
######################

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


### Define auxiliary functions
##############################

# def find_indefinable_n(p_value_cutoff):
#     print("[INFO] Looking for smallest number of read pairs n that yields a significant P-value with the given threshold of " + str(p_value_cutoff) + ".")
#     for n in range(1, 1000):
#         if dclass.calculate_binomial_p_value(n, 0) < p_value_cutoff:
#             print("\t[INFO] Smallest n: " + str(n) + " read pairs")
#             return n


def get_read_pair_count_paramters(diachromatic_interaction_file, n_indef_cutoff, p_value_cutoff):

    # find the smallest n for
    n_indef_cutoff = dclass.find_indefinable_n(p_value_cutoff)

    n_indefinable_list = []
    n_directed_list = []
    n_undirected_list = []
    n_undirected_reference_list = []

    n_dict = {}  # dictionary that stores the numbers of interactions with n read pairs

    with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
        n_interaction_total = 0
        line = fp.readline()
        while line:
            n_interaction_total += 1
            fields = line.rstrip('\n').split('\t')
            counts = fields[8].split(":")
            n_simple = int(counts[0])
            n_twisted = int(counts[1])
            n_total = n_simple + n_twisted
            if n_total < n_indef_cutoff:
                n_indefinable_list.append(n_total)
            else:
                p_val = dclass.calculate_binomial_p_value(n_simple, n_twisted)
                if p_val <= p_value_cutoff:
                    n_directed_list.append(n_total)
                    if n_total in n_dict:
                        n_dict[n_total] += 1
                    else:
                        n_dict[n_total] = 1
                else:
                    n_undirected_list.append(n_total)
            if n_interaction_total % 10000 == 0:
                print("\t[INFO]", n_interaction_total, "interactions processed ...")

            line = fp.readline()
    fp.close()

    # now use n_dict to select a reference set of interactions
    with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
        n_interaction_total = 0
        line = fp.readline()
        while line:
            n_interaction_total += 1
            fields = line.rstrip('\n').split('\t')
            counts = fields[8].split(":")
            n_simple = int(counts[0])
            n_twisted = int(counts[1])
            n_total = n_simple + n_twisted
            if n_total >= n_indef_cutoff:
                p_val = dclass.calculate_binomial_p_value(n_simple, n_twisted)
                if p_val >= p_value_cutoff:
                    if n_total in n_dict and  0 < n_dict[n_total]:
                        n_undirected_reference_list.append(n_total)
                        n_dict[n_total] = n_dict[n_total] - 1

            if n_interaction_total % 10000 == 0:
                print("\t[INFO]", n_interaction_total, "interactions processed ...")

            line = fp.readline()

    fp.close()

    # analyze arrays with read pair numbers

    mean_indefinable = np.mean(n_indefinable_list)
    mean_directed = np.mean(n_directed_list)
    mean_undirected = np.mean(n_undirected_list)
    mean_undirected_reference = np.mean(n_undirected_reference_list)

    sd_indefinable = np.std(n_indefinable_list)
    sd_directed = np.std(n_directed_list)
    sd_undirected = np.std(n_undirected_list)
    sd_undirected_reference = np.std(n_undirected_reference_list)


    median_indefinable = np.median(n_indefinable_list)
    median_directed = np.median(n_directed_list)
    median_undirected = np.median(n_undirected_list)
    median_undirected_reference = np.median(n_undirected_reference_list)

    print("\n")
    print("Mean read pair number - indefinable: " + str(mean_indefinable))
    print("Mean read pair number - directed: " + str(mean_directed))
    print("Mean read pair number - undirected: " + str(mean_undirected))
    print("Mean read pair number - undirected reference: " + str(mean_undirected_reference))
    print("\n")
    print("Standard deviation read pair number - indefinable: " + str(sd_indefinable))
    print("Standard deviation read pair number - directed: " + str(sd_directed))
    print("Standard deviation read pair number - undirected: " + str(sd_undirected))
    print("Standard deviation read pair number - undirected reference: " + str(sd_undirected_reference))
    print("\n")
    print("Median read pair number - indefinable: " + str(median_indefinable))
    print("Median read pair number - directed: " + str(median_directed))
    print("Median read pair number - undirected: " + str(median_undirected))
    print("Median read pair number - undirected reference: " + str(median_undirected_reference))

    # Create distance histograms
    plt.rcParams['figure.figsize'] = [20, 10]
    num_bins = 1000
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
    f.canvas.set_window_title('Distances between interacting digests')
    plt.title(out_prefix)
    plt.xlabel("Read pair numbers")
    plt.ylabel("Frequency")
    #plt.ylim(0, 600)

    ax1.set_title("Directed")
    n, bins, patches = ax1.hist(n_directed_list, num_bins, facecolor='blue', alpha=0.5, range=(0, 250))

    ax2.set_title("Undirected")
    n, bins, patches = ax2.hist(n_undirected_list, num_bins, facecolor='blue', alpha=0.5, range=(0, 250))

    ax3.set_title("Undirected reference")
    n, bins, patches = ax3.hist(n_undirected_reference_list, num_bins, facecolor='blue', alpha=0.5, range=(0, 250))

    ax4.set_title("Indefinable")
    n, bins, patches = ax4.hist(n_indefinable_list, num_bins, facecolor='blue', alpha=0.5, range=(0, 250))

    # Save the figure
    f.suptitle(out_prefix)
    figure_name = out_prefix + "_histogram.pdf"
    f.savefig(figure_name, bbox_inches='tight', format="pdf")


### Start execution
###################

#get_read_pair_count_paramters(diachromatic_interaction_file, 4, p_value_cutoff)

#exit(1)

# Init arrays for distribution of distances between digests
distance_array_simple = []
distance_array_twisted = []
distance_array_undirected = []
distance_array_undirected_reference = []
distance_array_indefinable = []

read_pair_num_array_simple = []
read_pair_num_array_twisted = []
read_pair_num_array_undirected = []
read_pair_num_array_undirected_reference = []
read_pair_num_array_indefinable = []

n_simple_interaction = 0
n_twisted_interaction = 0
n_undirected_interaction = 0
n_undirected_interaction_reference = 0
n_indefinable_interaction = 0

n_simple_interaction_rp = 0
n_twisted_interaction_rp = 0
n_undirected_interaction_rp = 0
n_undirected_interaction_reference_rp = 0
n_indefinable_interaction_rp = 0

n_trans_short_range_interaction = 0
n_non_status_pair_flag_interaction = 0

# Given the P-value cutoff, find the smallest n that yields a significant P-value
n_indefinable_cutoff = dclass.find_indefinable_n(p_value_cutoff)

# Determine distribution of n for directed interactions
min_digest_dist = 20000
n_dict = dclass.get_n_dict(diachromatic_interaction_file, status_pair_flag, min_digest_dist, p_value_cutoff)
n_dict_dir = n_dict.copy()
n_dict_ref = {}

# Iterate interactions and collect distances
print("[INFO] Determining distances between digest pairs in " + diachromatic_interaction_file + " ...")
with gzip.open(diachromatic_interaction_file, mode='rt') as fp:
    line = fp.readline()
    n_interaction_total = 0
    while line:

        n_interaction_total += 1

        # parse line representing one interaction
        interaction = dclass.Interaction(line)

        # restrict analysis to interactions between targeted promoters
        if status_pair_flag != "ALL" and interaction.get_digest_status_pair_flag() != status_pair_flag:
            n_non_status_pair_flag_interaction +=1
            line = fp.readline()
            continue

        # restrict analysis to cis long range interactions
        if not(interaction.is_cis_long_range(min_digest_dist)):
            n_trans_short_range_interaction += 1
            line = fp.readline()
            continue

        n_total = interaction.n_simple + interaction.n_twisted

        # Determine distance between the centers of digests
        distance = interaction.get_digest_distance()

        # set the type of interaction based on P-value ('S', 'T', 'U', 'NA')
        if interaction.get_interaction_type() == "TBD":
            interaction.set_interaction_type("TBD", p_value_cutoff, n_indefinable_cutoff)

        if interaction.get_interaction_type() == None:
            raise Exception("[FATAL] Interaction type is 'None'. This should never happen.")
        elif interaction.get_interaction_type() == "NA":
            distance_array_indefinable.append(distance)
            read_pair_num_array_indefinable.append(n_total)
            n_indefinable_interaction += 1
            n_indefinable_interaction_rp = n_indefinable_interaction_rp + n_total
        elif interaction.get_interaction_type() == "U":
            distance_array_undirected.append(distance)
            read_pair_num_array_undirected.append(n_total)
            n_undirected_interaction += 1
            n_undirected_interaction_rp = n_undirected_interaction_rp + n_total
            if n_total in n_dict and 0 < n_dict[n_total]:
                distance_array_undirected_reference.append(distance)
                read_pair_num_array_undirected_reference.append(n_total)
                n_undirected_interaction_reference += 1
                n_undirected_interaction_reference_rp = n_undirected_interaction_reference_rp + n_total
                n_dict[n_total] = n_dict[n_total] - 1
                if n_total in n_dict_ref:
                    n_dict_ref[n_total] += 1
                else:
                    n_dict_ref[n_total] = 1
        elif interaction.get_interaction_type() == "S":
            distance_array_simple.append(distance)
            read_pair_num_array_simple.append(n_total)
            n_simple_interaction += 1
            n_simple_interaction_rp = n_simple_interaction_rp + n_total
        elif interaction.get_interaction_type() == "T":
            distance_array_twisted.append(distance)
            read_pair_num_array_twisted.append(n_total)
            n_twisted_interaction += 1
            n_twisted_interaction_rp = n_twisted_interaction_rp + n_total
        else:
            line = fp.readline()
            print(interaction.get_interaction_type())
            raise Exception("[FATAL] Invalid interaction type. Should be either 'S', 'T', 'U' or 'NA' but was " + interaction.get_interaction_type() + ".")

        if n_interaction_total%100000 == 0:
            print("\t[INFO]", n_interaction_total, "interactions processed ...")

        line = fp.readline()

fp.close()

for key, value in n_dict.items():
    if 0 < value:
        print("[Warning] " + "Could not find corresponding number of undirected reference interactions for n = " + str(key) + ". Missing reference interactions: " + str(value))

# find largest n with significant interaction
n_max = -1
for i in n_dict:
    if n_max < i:
        n_max = i
print("Largest n with significant interaction: " + str(n_max))

print("[INFO] " + "Printing numbers of significant interactions with n read pairs to text file ...")
file_name = out_prefix + "_significant_and_reference_interactions.tab"
f_output = open(file_name, 'wt')
for n in range(0, n_max+1, 1):
    if n in n_dict:
        if n in n_dict_ref:
            print("n: " + str(n) + "\t" + str(n_dict_dir[n]) + "\t" + str(n_dict_ref[n]))
            f_output.write(str(n) + "\t" + str(n_dict_dir[n]) + "\t" + str(n_dict_ref[n]) + "\n")
        else:
            print("n: " + str(n) + "\t" + str(n_dict_dir[n]) + "\t" + "0")
            f_output.write(str(n) + "\t" + str(n_dict_dir[n]) + "\t" + "0")
    else:
        print("n: " + str(n) + "\t0\t0")
        f_output.write(str(n) + "\t0\t0")

f_output.close()

# Output summary
print("[INFO] " + "Summary statistics")
print("\t[INFO] Total number of interactions: " + str(n_interaction_total))
print("\t[INFO] Number of non " + status_pair_flag + " interactions: " + str(n_non_status_pair_flag_interaction) + " (discarded)")
print("\t[INFO] Number of trans and short range interactions: " + str(n_trans_short_range_interaction) + " (discarded)")
print("\t------------------------------------------------------------------")
print("\t[INFO] Number of simple interactions: " + str(n_simple_interaction))
print("\t[INFO] Number of twisted interactions: " + str(n_twisted_interaction))
print("\t[INFO] Number of undirected interactions: " + str(n_undirected_interaction))
print("\t[INFO] Number of undirected interactions reference: " + str(n_undirected_interaction_reference))
print("\t[INFO] Number of indefinable interactions: " + str(n_indefinable_interaction))
print("\t------------------------------------------------------------------")
print("\t[INFO] Number of read pairs in simple interactions: " + str(n_simple_interaction_rp) + " (" + str(n_simple_interaction_rp/n_simple_interaction) + ")")
print("\t[INFO] Number of read pairs in twisted interactions: " + str(n_twisted_interaction_rp) + " (" + str(n_twisted_interaction_rp/n_twisted_interaction) + ")")
print("\t[INFO] Number of read pairs in undirected interactions: " + str(n_undirected_interaction_rp) + " (" + str(n_undirected_interaction_rp/n_undirected_interaction) + ")")
print("\t[INFO] Number of read pairs in undirected interactions reference: " + str(n_undirected_interaction_reference_rp) + " (" + str(n_undirected_interaction_reference_rp/n_undirected_interaction_reference) + ")")
print("\t[INFO] Number of read pairs in indefinable interactions: " + str(n_indefinable_interaction_rp) + " (" + str(n_indefinable_interaction_rp/n_indefinable_interaction) + ")")
print("\t------------------------------------------------------------------")
# print("[INFO] " + "Writing numpy arrays with distances to disk ...")
# file_path_name = out_prefix + "_distance_array_simple"
# np.save(file_path_name, np.array(distance_array_simple))
# file_path_name = out_prefix + "_distance_array_twisted"
# np.save(file_path_name, np.array(distance_array_twisted))
# file_path_name = out_prefix + "_distance_array_undirected"
# np.save(file_path_name, np.array(distance_array_undirected))
# file_path_name = out_prefix + "_distance_array_undirected_reference"
# np.save(file_path_name, np.array(distance_array_undirected_reference))
# file_path_name = out_prefix + "_distance_array_indefinable"
# np.save(file_path_name , np.array(distance_array_indefinable))

print("[INFO] " + "Writing arrays with distances to files ...")

file_name = out_prefix + "_digest_distances_simple.txt"
f_output = open(file_name, 'wt')
for i in range(len(distance_array_simple)):
   f_output.write(str(distance_array_simple[i]) + "\t")
   f_output.write(str(read_pair_num_array_simple[i]) + "\n")
f_output.close()

file_name = out_prefix + "_digest_distances_twisted.txt"
f_output = open(file_name, 'wt')
for i in range(len(distance_array_twisted)):
   f_output.write(str(distance_array_twisted[i]) + "\t")
   f_output.write(str(read_pair_num_array_twisted[i]) + "\n")
f_output.close()

file_name = out_prefix + "_digest_distances_undirected.txt"
f_output = open(file_name, 'wt')
for i in range(len(distance_array_undirected)):
   f_output.write(str(distance_array_undirected[i]) + "\t")
   f_output.write(str(read_pair_num_array_undirected[i]) + "\n")
f_output.close()

file_name = out_prefix + "_digest_distances_undirected_reference.txt"
f_output = open(file_name, 'wt')
for i in range(len(distance_array_undirected_reference)):
   f_output.write(str(distance_array_undirected_reference[i]) + "\t")
   f_output.write(str(read_pair_num_array_undirected_reference[i]) + "\n")
f_output.close()

file_name = out_prefix + "_digest_distances_indefinable.txt"
f_output = open(file_name, 'wt')
for i in range(len(distance_array_indefinable)):
   f_output.write(str(distance_array_indefinable[i]) + "\t")
   f_output.write(str(read_pair_num_array_indefinable[i]) + "\n")
f_output.close()

print("[INFO] " + "... done.")
