#!/usr/bin/env python
import diachrscripts_classes as dclass
import argparse
import gzip
import statistics as statistics
import matplotlib.pyplot as plt

# Parse command line arguments
parser = argparse.ArgumentParser(description='Analyze distances between interacting digest pairs.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--interaction-file', help='Diachromatic interaction file.')
parser.add_argument('--status-pair-flag', help='Pair of \'A\' and \'I\' depending on whether a digest was selected for enrichment (A) or not (I).')

args = parser.parse_args()
out_prefix = args.out_prefix
diachromatic_interaction_file = args.interaction_file
status_pair_flag = args.status_pair_flag
if status_pair_flag != "ALL":
    status_pair_flag = sorted(status_pair_flag)[0] + sorted(status_pair_flag)[1]

# Init arrays for distribution of distances between digests
distance_array_simple = []
distance_array_twisted = []
distance_array_undirected = []

n_significant_simple_interactions = 0
n_significant_twisted_interactions = 0
n_undirected_interactions = 0

# Iterate interactions and collect distances
with gzip.open(diachromatic_interaction_file, 'r' + 't') as fp:
    line = fp.readline()
    while line:

        # Parse line
        values = line.split("\t")

        digest_1 = dclass.Digest(values[0], int(values[1]), int(values[2]))
        if values[3] == 'A':
            digest_1.set_active()
        digest_2 = dclass.Digest(values[4], int(values[5]), int(values[6]))
        if values[7] == 'A':
            digest_2.set_active()

        values2 = values[8].split(":")
        n_simple = int(values2[0])
        n_twisted = int(values2[1])
        interaction = dclass.Interaction(digest_1, digest_2, n_simple, n_twisted)

        # Restrict analysis to subset of interactions, e.g. 'AA'
        if(status_pair_flag != "ALL" and status_pair_flag != interaction.get_status_pair_flag()):
            line = fp.readline()
            continue

        # Restrict analysis to cis interactions
        if not(interaction.is_cis()):
            line = fp.readline()
            continue

        # Get binomial P-value
        p_value = interaction.get_binomial_p_value()

        # Determine distance between the centers of digests
        distance = interaction.get_digest_distance()

        # Collect distances
        if p_value <= 0.05: # significant interactions

            # simple
            if n_twisted < n_simple:
                distance_array_simple.append(distance)
                n_significant_simple_interactions += 1
            else: # twisted
                distance_array_twisted.append(distance)
                n_significant_twisted_interactions += 1

        else: # undirected
            distance_array_undirected.append(distance)
            n_undirected_interactions += 1

        line = fp.readline()

fp.close()

# Output summary
print("Summary")
print("=======\n")
print("Analysis for: " + out_prefix)
print("Input file: " + diachromatic_interaction_file)
print("Status pair flag: " + status_pair_flag)
print("\n\n")

print("Number of significant simple interactions: " + str(n_significant_simple_interactions))
print("Number of significant twisted simple interactions: " + str(n_significant_twisted_interactions))
print("Number of undirected interactions: " + str(n_significant_twisted_interactions))
print("\n\n")

# Determine mean distances
mean_simple = statistics.mean(distance_array_simple)
print("Mean simple: " + str(mean_simple))
mean_twisted = statistics.mean(distance_array_twisted)
print("Mean twisted: " + str(mean_twisted))
mean_undirected = statistics.mean(distance_array_undirected)
print("Mean undirected: " + str(mean_undirected))
print("\n")

# Determine median distances
median_simple = statistics.median(distance_array_simple)
print("Median simple: " + str(median_simple))
median_twisted = statistics.median(distance_array_twisted)
print("Median twisted: " + str(median_twisted))
median_undirected = statistics.median(distance_array_undirected)
print("Median undirected: " + str(median_undirected))

# Create plot
num_bins = 100
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)
f.canvas.set_window_title('Distances between interacting digests')
plt.xlabel("Distance")
plt.ylabel("Frequency")

ax1.set_title("Simple")
n, bins, patches = ax1.hist(distance_array_simple, num_bins, facecolor='blue', alpha=0.5, range=(0,100000))

ax2.set_title("Twisted")
n, bins, patches = ax2.hist(distance_array_twisted, num_bins, facecolor='blue', alpha=0.5, range=(0,100000))

ax3.set_title("Undirected")
n, bins, patches = ax3.hist(distance_array_undirected, num_bins, facecolor='blue', alpha=0.5, range=(0,100000))

plt.show()
