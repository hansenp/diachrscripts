#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import statistics as statistics
from scipy.stats import binom
import gzip

# Arguments
interaction_file = sys.argv[1] # Diachromatic interaction file gzipped
interaction_category = sys.argv[2] # Interaction category (TBA)

distance_array_simple = []
distance_array_twisted = []
distance_array_undirected = []

with gzip.open(interaction_file, 'r' + 't') as fp:
    line = fp.readline()
    while line:
        values = line.split("\t")

        chr_name_1 = values[0]
        d_sta_1 = int(values[1])
        d_end_1 = int(values[2])
        d_state_1 = values[3]

        chr_name_2 = values[4]
        d_sta_2 = int(values[5])
        d_end_2 = int(values[6])
        d_state_2 = values[7]

        values2 = values[8].split(":")
        simple_cnt = int(values2[0])
        twisted_cnt = int(values2[1])

        # restrict analysis to subset of interactions
        if (interaction_category != d_state_1 + d_state_2) & (interaction_category != d_state_2 + d_state_1) & (
                interaction_category != "ALL"):
            line = fp.readline()
            continue

        # restrict analysis to cis interactions
        if chr_name_1 != chr_name_2:
            line = fp.readline()
            continue

        # calculate binomial P-value
        if simple_cnt < twisted_cnt:
            p_value = 1 - binom.cdf(twisted_cnt-1, simple_cnt + twisted_cnt, 0.5)
        else:
            p_value = 1 - binom.cdf(simple_cnt-1, simple_cnt + twisted_cnt, 0.5)

        # determine distance between the centers of digests
        center_1 = int(d_sta_1 + (d_end_1 - d_sta_1) / 2)
        center_2 = int(d_sta_2 + (d_end_2 - d_sta_2) / 2)
        distance = center_2 - center_1

        if p_value <= 0.05: # significant interactions

            # simple
            if twisted_cnt < simple_cnt:
                distance_array_simple.append(distance)
            else: # twisted
                distance_array_twisted.append(distance)

        else: #undirected
            distance_array_undirected.append(distance)

        line = fp.readline()

fp.close()

# determine mean distances
mean_simple = statistics.mean(distance_array_simple)
print "Mean simple:", mean_simple
mean_twisted = statistics.mean(distance_array_twisted)
print "Mean twisted:", mean_twisted
mean_undirected = statistics.mean(distance_array_undirected)
print "Mean undirected:", mean_undirected

# determine median distances
median_simple = statistics.median(distance_array_simple)
print "Median simple:", median_simple
median_twisted = statistics.median(distance_array_twisted)
print "Median twisted:", median_twisted
median_undirected = statistics.median(distance_array_undirected)
print "Median undirected:", median_undirected

# create plot
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



