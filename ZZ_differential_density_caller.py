"""
This script takes as input two BedGraph files with digest densities created with the script
'diachrscripts/ZZ_get_digest_density_profile_bedgraphs.py', calculates the local correlation within a sliding
window with a specified number of digests (default is 25), and creates a new BedGraph file with the local correlations.
"""

import argparse
import numpy
import scipy.stats

### Parse command line
######################

parser = argparse.ArgumentParser(description='Determine differential digest densities along the genome.')
parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--bed-graph-file-1', help='BedGraph file with digest densities.', required=True)
parser.add_argument('--bed-graph-file-2', help='Another BedGraph file with digest densities.', required=True)
parser.add_argument('--bed-graph-file-bl', help='Another BedGraph file with densities of blacklisted digests.', required=True)
parser.add_argument('--data-set-track-tag', help='A tag that is used to create track names and descriptions.')
parser.add_argument('--digest-density-range', help='Number of consecutive digests in the sliding window.', default=25)
parser.add_argument('--correlation-p-value-cutoff', help='Digests with smaller P-values will be written to a BED file.', default=0.0001)

args = parser.parse_args()
out_prefix = args.out_prefix
bed_graph_file_1 = args.bed_graph_file_1
bed_graph_file_2 = args.bed_graph_file_2
bed_graph_file_bl = args.bed_graph_file_bl
data_set_track_tag = args.data_set_track_tag
digest_density_range = int(args.digest_density_range)
correlation_p_value_cutoff = float(args.correlation_p_value_cutoff)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --bed-graph-file-1: " + bed_graph_file_1)
print("\t[INFO] --bed-graph-file-2: " + bed_graph_file_2)
print("\t[INFO] --bed-graph-file-bl: " + bed_graph_file_bl)
print("\t[INFO] --data-set-track-tag: " + data_set_track_tag)
print("\t[INFO] --digest-density-range: " + str(digest_density_range))
print("\t[INFO] --correlation-p-value-cutoff: " + '{:10.5f}'.format(correlation_p_value_cutoff))

# Define UCSC track options
track_maxHeightPixels = 50

track_color = "255,161,0"
track_altColor ="0,100,200"
track_name = data_set_track_tag + " - DIGEST DENSITY DIFFERENCE"
track_description = track_name
bedgraph_stream_digests_density_difference_output = open(out_prefix + "_digests_density_difference.bedgraph", 'wt')
bedgraph_stream_digests_density_difference_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_color = "255,0,0"
track_altColor ="0,0,255"
track_name = data_set_track_tag + " - DIGEST DENSITY CORRELATION - r"
track_description = track_name
bedgraph_stream_digests_density_correlation_r_output = open(out_prefix + "_digests_density_correlation_r.bedgraph", 'wt')
bedgraph_stream_digests_density_correlation_r_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_name = data_set_track_tag + " - DIGEST DENSITY CORRELATION - P-val"
track_description = track_name
bedgraph_stream_digests_density_correlation_p_output = open(out_prefix + "_digests_density_correlation_p.bedgraph", 'wt')
bedgraph_stream_digests_density_correlation_p_output.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')


track_color = "70,125,70"
track_name = data_set_track_tag + " - DIGEST DENSITY CORRELATION - SIG 1"
track_description = track_name
bedgraph_stream_digests_density_correlation_sig1_output = open(out_prefix + "_digests_density_correlation_sig1.bed", 'wt')
bedgraph_stream_digests_density_correlation_sig1_output.write("track type=bed name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')

track_color = "70,70,125  "
track_name = data_set_track_tag + " - DIGEST DENSITY CORRELATION - SIG 2"
track_description = track_name
bedgraph_stream_digests_density_correlation_sig2_output = open(out_prefix + "_digests_density_correlation_sig2.bed", 'wt')
bedgraph_stream_digests_density_correlation_sig2_output.write("track type=bed name=\"" + track_name + "\" description=\"" + track_description + "\" visibility=full color=" + track_color + " altColor=" + track_altColor + " priority=20 maxHeightPixels=" + str(track_maxHeightPixels) + '\n')


### Iterate BedGrap files
#########################



bg_file_1  = open(bed_graph_file_1, 'rt')
bg_file_2  = open(bed_graph_file_2, 'rt')
bg_file_bl  = open(bed_graph_file_bl, 'rt')

bg1_line = bg_file_1.readline().strip('\n')
bg1_line = bg_file_1.readline().strip('\n') # Skip header line
bg2_line = bg_file_2.readline().strip('\n')
bg2_line = bg_file_2.readline().strip('\n') # Skip header line
bgbl_line = bg_file_bl.readline().strip('\n')
bgbl_line = bg_file_bl.readline().strip('\n') # Skip header line

n_progress = 0

digest_queue_1 = []
digest_queue_2 = []
digest_queue_bl = []
r = 0.0
p = 0.0
chr_d_cnt = 0

last_chr = "chr1"

while bg1_line and bg2_line:

    # Report progress
    n_progress += 1
    if n_progress % 100000 == 0:
        print("\tProcessed " + str(n_progress) + " digests ...")

    (chr_1, sta_1, end_1, dens_1) = bg1_line.split('\t')
    (chr_2, sta_2, end_2, dens_2) = bg2_line.split('\t')
    (chr_2, sta_2, end_2, dens_bl) = bgbl_line.split('\t')

    dens_1 = float(dens_1)
    dens_2 = float(dens_2)
    dens_bl = float(dens_bl)

    # Set everything to zero when encountering a new chromosome
    if last_chr != chr_1:

        # Empty queue
        while 0 < len(digest_queue_1):
            digest_queue_1.pop(0)
            digest_queue_2.pop(0)
        chr_d_cnt = 0

        # Reset variables for sliding digest window
        last_chr = chr_1
        chr_d_cnt = 0
        r = 0.0

    # Fill queue
    digest_queue_1.append(dens_1)
    digest_queue_2.append(dens_2)
    chr_d_cnt += 1

    # Update digest densities
    if len(digest_queue_1) == digest_density_range + 1:

        # Remove old digest from queue
        digest_queue_1.pop(0)
        digest_queue_2.pop(0)

        # Calculate correlation
        r, p = scipy.stats.pearsonr(numpy.array(digest_queue_1).astype(numpy.float), numpy.array(digest_queue_2).astype(numpy.float))

    if not numpy.isnan(r):
        x = r
        y = numpy.emath.log10(p)
    else:
        x = 2.0
        y = 2.0

        # Keep track of current chromosome
        last_chr = chr_1

    if y < correlation_p_value_cutoff and dens_2 < dens_1:
        bedgraph_stream_digests_density_correlation_sig1_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(y) + ':' + str(dens_1) + ':' + str(dens_2) + ':' + str(r) + '\n')

    if y < correlation_p_value_cutoff and dens_2 >= dens_1:
        bedgraph_stream_digests_density_correlation_sig2_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(y) + ':' + str(dens_1) + ':' + str(dens_2) + ':' + str(r) + '\n')

    # Center window
    window_center_pos = int(digest_density_range / 2) - 1 + digest_density_range % 2  # works for uneven and even window sizes

    # Write window centered correlations and P-values to file
    if chr_d_cnt - 1 < window_center_pos:
        bedgraph_stream_digests_density_correlation_r_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(0.0) + '\n')
        bedgraph_stream_digests_density_correlation_p_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(0.0) + '\n')
    elif digest_density_range <= chr_d_cnt:
        bedgraph_stream_digests_density_correlation_r_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(x) + '\n')
        bedgraph_stream_digests_density_correlation_p_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(y) + '\n')

    if 0 < dens_2:
        bedgraph_stream_digests_density_difference_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(1.0-dens_1/(dens_2+dens_bl)    ) + '\n')
    else:
        bedgraph_stream_digests_density_difference_output.write(chr_1 + '\t' + sta_1 + '\t' + end_1 + '\t' + str(-1.0) + '\n')

    bg1_line = bg_file_1.readline().strip('\n')
    bg2_line = bg_file_2.readline().strip('\n')


bg_file_1.close()
bg_file_2.close()

bedgraph_stream_digests_density_difference_output.close()
bedgraph_stream_digests_density_correlation_r_output.close()
bedgraph_stream_digests_density_correlation_p_output.close()
bedgraph_stream_digests_density_correlation_sig1_output.close()
bedgraph_stream_digests_density_correlation_sig2_output.close()
