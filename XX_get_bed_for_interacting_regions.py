"""
This script takes as input an enhanced interaction file and determines two digest sets:

   1. Digests that are involved in at least one directed interaction.
   2. Digest that are involved in undirected reference interactions and in no directed interaction.

The enhanced interaction file is then iterated a second time and the digest regions of each interaction are written
to one of two BED files that:

   1. Contain digest regions that are involved in at least one directed interaction.
   2. Contain digest regions that are involved in undirected reference interactions and in no directed interaction.

Since one digest can be involved in more than one interaction, certain regions may occur in the BED file multiple times.
"""


import argparse
import gzip
import diachrscripts_toolkit as dclass
import numpy


### Parse command line
######################

parser = argparse.ArgumentParser(description='Get distances to nearest TAD boundary for diegstes involved in directed interactions and undirected reference interactions not involved in directed interactions.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enhanced-interaction-file', help='Enhanced interaction file supplemented with digest associated gene symbols and TSS as well as directionality P-values.', required=True)
parser.add_argument('--gopher-digest-file', help='File that contains the coordinates of all digests for a given genome build.', required=True)

args = parser.parse_args()
out_prefix = args.out_prefix
enhanced_interaction_file = args.enhanced_interaction_file
gopher_digest_file = args.gopher_digest_file

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enhanced-interaction-file: " + enhanced_interaction_file)
print("\t[INFO] --gopher-digest-file: " + gopher_digest_file)

# Init digest sets
digests_from_dir_inter = set()
digests_from_undir_inter = set()

print("[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        # Add digest to set for directed interactions
        if interaction_category == "DIII":
            digests_from_dir_inter.add(chr_a + '\t' + str(sta_a) + '\t' + str(end_a))
            digests_from_dir_inter.add(chr_b + '\t' + str(sta_b) + '\t' + str(end_b))

        # Add digest to set for undirected interactions
        if interaction_category == "UIRII":
            digests_from_undir_inter.add(chr_a + '\t' + str(sta_a) + '\t' + str(end_a))
            digests_from_undir_inter.add(chr_b + '\t' + str(sta_b) + '\t' + str(end_b))

        line = fp.readline()

fp.close()

digests_from_undir_inter_wo_digests_from_dir_inter = digests_from_undir_inter - digests_from_dir_inter
digests_from_undir_inter_wo_digests_from_dir_inter = digests_from_undir_inter

print("\tNumber of digests from directed interactions: " + str(len(digests_from_dir_inter)))
print("\tNumber of digests from undirected interactions: " + str(len(digests_from_undir_inter_wo_digests_from_dir_inter)))

print("[INFO] ... done.")

directed_interaction_ends = dict()
undirected_interaction_ends = dict()

directed_interaction_ends_distances = dict()
undirected_interaction_ends_distances = dict()

########################################
# Get BED file with interacting digests

bed_stream_dir_inter_digests_output = open(out_prefix + "_dir_inter_digests.bed", 'wt')
bed_stream_undir_inter_digests_output = open(out_prefix + "_undir_inter_digests.bed", 'wt')

n_directed_interaction_ends = 0
n_undirected_interaction_ends = 0

print("[INFO] Iterating enhanced interaction file ...")
with gzip.open(enhanced_interaction_file, 'rt') as fp:

    n_interaction_total = 0
    line = fp.readline()
    while line:

        # Parse enhanced interactions line
        chr_a, sta_a, end_a, syms_a, tsss_a, chr_b, sta_b, end_b, syms_b, tsss_b, enrichment_pair_tag, strand_pair_tag, interaction_category, neg_log_p_value, rp_total, i_dist = \
            dclass.parse_enhanced_interaction_line_with_gene_symbols(line)

        x = line.split('\t')
        n_simple = int(x[4].split(':')[0])
        n_twisted = int(x[4].split(':')[1])
        n = n_simple + n_twisted

        if interaction_category == "DIII":

            if n_simple < n_twisted:
                bed_name = str(n_simple) + ':' + str(n_twisted) + ';' +  str(neg_log_p_value) + ';' + str(i_dist) + ";T"
            else:
                bed_name = str(n_simple) + ':' + str(n_twisted) + ';' +  str(neg_log_p_value) + ';' + str(i_dist) + ";S"

            # Write interaction end to BED file
            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) in digests_from_dir_inter:
                bed_name_b = chr_b + ':' + str(sta_b) + '-' + str(end_b) + ';' + bed_name
                bed_stream_dir_inter_digests_output.write(chr_a + '\t' + str(sta_a) + '\t' + str(end_a) + '\t' + bed_name_b + '\n')
                n_directed_interaction_ends += 1
                if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) not in directed_interaction_ends:
                    directed_interaction_ends[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] = 1
                    for i in range(1,n):
                        directed_interaction_ends_distances[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] = [i_dist]
                else:
                    directed_interaction_ends[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] += 1
                    for i in range(1,n):
                        directed_interaction_ends_distances[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)].append(i_dist)

            # Write interaction end to BED file
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) in digests_from_dir_inter:
                bed_name_a = chr_a + ':' + str(sta_a) + '-' + str(end_a) + ';' + bed_name
                bed_stream_dir_inter_digests_output.write(chr_b + '\t' + str(sta_b) + '\t' + str(end_b) + '\t' + bed_name_a + '\n')
                n_directed_interaction_ends += 1
                if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) not in directed_interaction_ends:
                    directed_interaction_ends[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] = 1
                    for i in range(1,n):
                        directed_interaction_ends_distances[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] = [i_dist]
                else:
                    directed_interaction_ends[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] += 1
                    for i in range(1,n):
                        directed_interaction_ends_distances[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)].append(i_dist)

        if interaction_category == "UIRII":

            bed_name = str(n_simple) + ':' + str(n_twisted) + ';' +  str(neg_log_p_value) + ';' + str(i_dist) + ";U"

            # Add digest to set for directed interactions
            if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) in digests_from_undir_inter_wo_digests_from_dir_inter:
                bed_name_b = chr_b + ':' + str(sta_b) + '-' + str(end_b) + ';' + bed_name
                bed_stream_undir_inter_digests_output.write(chr_a + '\t' + str(sta_a) + '\t' + str(end_a) + '\t' + bed_name_b + '\n')
                n_undirected_interaction_ends += 1
                if chr_a + '\t' + str(sta_a) + '\t' + str(end_a) not in undirected_interaction_ends:
                    undirected_interaction_ends[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] = 1
                    for i in range(1,n):
                        undirected_interaction_ends_distances[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] = [i_dist]
                else:
                    undirected_interaction_ends[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)] += 1
                    for i in range(1,n):
                        undirected_interaction_ends_distances[chr_a + '\t' + str(sta_a) + '\t' + str(end_a)].append(i_dist)

            # Add digest to set for undirected interactions
            if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) in digests_from_undir_inter_wo_digests_from_dir_inter:
                bed_name_a = chr_a + ':' + str(sta_a) + '-' + str(end_a) + ';' + bed_name
                bed_stream_undir_inter_digests_output.write(chr_b + '\t' + str(sta_b) + '\t' + str(end_b) + '\t' + bed_name_a + '\n')
                n_undirected_interaction_ends += 1
                if chr_b + '\t' + str(sta_b) + '\t' + str(end_b) not in undirected_interaction_ends:
                    undirected_interaction_ends[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] = 1
                    for i in range(1,n):
                        undirected_interaction_ends_distances[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] = [i_dist]
                else:
                    undirected_interaction_ends[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)] += 1
                    for i in range(1,n):
                        undirected_interaction_ends_distances[chr_b + '\t' + str(sta_b) + '\t' + str(end_b)].append(i_dist)

        line = fp.readline()

fp.close()
bed_stream_dir_inter_digests_output.close()
bed_stream_undir_inter_digests_output.close()


print("\tNumber of directed interactions ending in directed digests: " + str(n_directed_interaction_ends))
print("\tNumber of undirected interactions ending in undirected digests: " + str(n_undirected_interaction_ends))

print("\tPile-up coefficient directed: " + str(len(digests_from_dir_inter)/n_directed_interaction_ends))
print("\tPile-up coefficient undirected: " + str(len(digests_from_undir_inter_wo_digests_from_dir_inter)/n_undirected_interaction_ends))

print()
print("out_prefix\tdig_from_dir_inter\tdig_from_undir_inter\tdig_from_undir_inter_wo_dig_from_dir_inter\tn_directed_interaction_ends\tn_undirected_interaction_ends")
print(out_prefix + "\t" + str(len(digests_from_dir_inter)) + "\t" + str(len(digests_from_undir_inter)) + '\t' + str(len(digests_from_undir_inter_wo_digests_from_dir_inter)) + '\t' + str(n_directed_interaction_ends) + '\t' + str(n_undirected_interaction_ends))

bedgraph_stream_dir_inter_digests_output = open(out_prefix + "_dir_inter_digests.bedgraph", 'wt')
bedgraph_stream_dir_inter_digests_output.write("track type=bedGraph name=\"Directed\" description=\"Directed\" visibility=full color=255,161,0 altColor=0,100,200 priority=20" + '\n')

bedgraph_stream_undir_inter_digests_output = open(out_prefix + "_undir_inter_digests.bedgraph", 'wt')
bedgraph_stream_undir_inter_digests_output.write("track type=bedGraph name=\"Undirected\" description=\"Undirected\" visibility=full color=071,115,130 altColor=0,100,200 priority=20" + '\n')

bedgraph_avg_dist_stream_dir_inter_digests_output = open(out_prefix + "_dir_inter_digests_i_dist.bedgraph", 'wt')
bedgraph_avg_dist_stream_dir_inter_digests_output.write("track type=bedGraph name=\"Directed - dist -rp\" description=\"Directed - dist - rp\" visibility=full color=255,161,0 altColor=0,100,200 priority=20" + '\n')

bedgraph_avg_dist_stream_undir_inter_digests_output = open(out_prefix + "_undir_inter_digests_i_dist.bedgraph", 'wt')
bedgraph_avg_dist_stream_undir_inter_digests_output.write("track type=bedGraph name=\"Undirected - dist - rp\" description=\"Undirected - dist - rp\" visibility=full color=0,100,200 altColor=0,100,200 priority=20" + '\n')


print("[INFO] Iterating digest file ...")
with open(gopher_digest_file, 'rt') as fp:
    line = fp.readline()
    line = fp.readline()
    while line:
        fields = line.split('\t')
        coord_key = fields[0] + '\t' + fields[1] + '\t' + fields[2]
        if coord_key in directed_interaction_ends:
            dir_inter_ends = directed_interaction_ends[coord_key]
        else:
            dir_inter_ends = 0

        if coord_key in undirected_interaction_ends:
            undir_inter_ends = undirected_interaction_ends[coord_key]
        else:
            undir_inter_ends = 0

        # calculate score and write to BedGraph for directed
        if 0 < dir_inter_ends:
            directionality_score = dir_inter_ends * (dir_inter_ends / (dir_inter_ends + undir_inter_ends))
            median_i_dist_dir = numpy.median(directed_interaction_ends_distances[coord_key])
        else:
            directionality_score = 0
            median_i_dist_dir = 0

        bedgraph_stream_dir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1])-1) + '\t' + str(fields[2]) + '\t' + str(directionality_score) + '\n')

        bedgraph_avg_dist_stream_dir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1])-1) + '\t' + str(fields[2]) + '\t' + str(median_i_dist_dir) + '\n')


        # calculate score and write to BedGraph for undirected
        if 0 < undir_inter_ends:
            directionality_score = undir_inter_ends * (undir_inter_ends / (undir_inter_ends + dir_inter_ends))
            median_i_dist_undir = numpy.median(undirected_interaction_ends_distances[coord_key])
        else:
            directionality_score = 0
            median_i_dist_undir = 0

        bedgraph_stream_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1])-1) + '\t' + str(fields[2]) + '\t' + str(directionality_score) + '\n')

        bedgraph_avg_dist_stream_undir_inter_digests_output.write(fields[0] + '\t' + str(int(fields[1])-1) + '\t' + str(fields[2]) + '\t' + str(median_i_dist_undir) + '\n')

        line = fp.readline()

fp.close()

bedgraph_stream_dir_inter_digests_output.close()
bedgraph_stream_undir_inter_digests_output.close()

bedgraph_avg_dist_stream_dir_inter_digests_output.close()
bedgraph_avg_dist_stream_undir_inter_digests_output.close()

exit(0)
for key in directed_interaction_ends:
    if key in undirected_interaction_ends:
        if undirected_interaction_ends[key] < directed_interaction_ends[key]:
            print(key)
            print(directed_interaction_ends[key])
            print(str(undirected_interaction_ends[key]))