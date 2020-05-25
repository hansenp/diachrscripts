#!/usr/bin/env python

"""
This script takes a tab separated file with the interactions analyzed in paper of Javierre et al. (2016) and creates
one BED file for each of the 17 human primary hematopoietic cell types with the start coordinates of the first and the
end coordinates of the second digest in sequential order. For each cell type, interactions with a CHiCAGO score
lower or a chosen threshold (default is 5) will be discarded.


The input file has the following structure:

* Field 1: baitChr - Chromosome of this baited digest.

* Field 2: baitStart - Start coordinate of this baited digest.

* Field 3: baitEnd - End coordinate of this baited digest.

* Field 4: baitID - Restriction fragment ID of this baited digest.

* Field 5: baitName - Names of genes whose promoters map to this baited digest.

* Field 6: oeChr - Chromosome of the other end (oe) digest that may also be baited.

* Field 7: oeStart - Start coordinate of the other end digest.

* Field 8: oeEnd - End coordinate of the other end digest.

* Field 9: oeID - Restriction fragment ID of the other end digest.

* Field 10: oeName - Restriction fragment ID of the other end digest. If this is not ‘.’ this interaction is a promoter-promoter interaction.

* Field 11: dist - Distance between baited digest and other end digest.
    - Negative, if other end digest comes before baited digest in sequential order
    - NA, if this is a trans-chromosomal interaction

* Fields 12 to 28: Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, tB

* Fields 29 and 30: Not relevant.


Each of the 17 output BED files has the following structure:

* Field 1: Chromosome of the first digest in sequential order.

* Field 2: Start coordinate of the first digest of an interaction in sequential order.

* Field 3: End coordinate of the second digest of an interaction in sequential order.


There will be one additional output file with the number of discarded and kept interactions for each cell type.

"""


import argparse
import gzip


### Parse command line
######################

parser = argparse.ArgumentParser(description='From the tsv file \'PCHiC_peak_matrix_cutoff5.tsv\', create one BED file for each of the 17 cell types that contains interactions that have a CHiCAGO equal or greater 5.')

parser.add_argument('--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--tsv-interaction-file', help='tsv file taken from the supplementary materials of the Javierre paper.', required=True)
parser.add_argument('--chicago-threshold', help='Use this switch to estimate a P-value cutoff that corresponds to a given FDR threshold.', default=5.0)

args = parser.parse_args()
out_prefix = args.out_prefix
tsv_interaction_file = args.tsv_interaction_file
chicago_threshold = float(args.chicago_threshold)

print("[INFO] " + "Input parameters")
print("\t[INFO] --out_prefix: " + out_prefix)
print("\t[INFO] --tsv-interaction-file: " + tsv_interaction_file)
print("\t[INFO] --chicago-threshold: " + str(chicago_threshold))


### Start execution
###################

print("[INFO] Iterating tsv file ...")
n_interaction = 0
with open(tsv_interaction_file, 'r' + 't') as fp:

    for line in fp:

        # Count total number of interactions
        n_interaction += 1

        # Report progress
        if n_interaction % 100000 == 0:
            print("\t\t[INFO]", n_interaction, "interactions processed ...")

        # Split line of tsv file
        baitChr, baitStart, baitEnd, baitID, baitName, oeChr, oeStart, oeEnd, oeID, oeName, dist, Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, tB, clusterID, clusterPostProb = line.split("\t")

        # If this is the first line, check whether everything is where it should be
        if(n_interaction == 1):
            if("baitChr" != baitChr or "baitStart" != baitStart or "baitEnd" != baitEnd or "baitID" != baitID or
                    "baitName" != baitName or "oeChr" != oeChr or "oeStart" != oeStart or "oeEnd" != oeEnd or
                    "oeID" != oeID or "oeName" != oeName or "dist" != dist or "Mon" != Mon or "Mac0" != Mac0 or
                    "Mac1" != Mac1 or "Mac2" != Mac2 or "Neu" != Neu or "MK" != MK or "EP" != EP or "Ery" != Ery or
                    "FoeT" != FoeT or "nCD4" != nCD4 or "tCD4" != tCD4 or "aCD4" != aCD4 or "naCD4" != naCD4 or
                    "nCD8" != nCD8 or "tCD8" != tCD8):
                print("baitChr: " + baitChr)
                print("baitStart: " + baitStart)
                print("baitEnd: " + baitEnd)
                print("baitID: " + baitID)
                print("baitName: " + baitName)
                print("oeChr: " + oeChr)
                print("oeStart: " + oeStart)
                print("oeEnd: " + oeEnd)
                print("oeID: " + oeID)
                print("oeName: " + oeName)
                print("dist: " + dist)
                print("Mon: " + Mon)
                print("Mac0: " + Mac0)
                print("Mac1: " + Mac1)
                print("Mac2: " + Mac2)
                print("Neu: " + Neu)
                print("MK: " + MK)
                print("EP: " + EP)
                print("Ery: " + Ery)
                print("FoeT: " + FoeT)
                print("nCD4: " + nCD4)
                print("tCD4: " + tCD4)
                print("aCD4: " + aCD4)
                print("naCD4: " + naCD4)
                print("nCD8: " + nCD8)
                print("tCD8: " + tCD8)
                print("nB: " + nB)
                print("tB: " + tB)
                print("clusterID: " + clusterID)
                print("clusterPostProb: " + clusterPostProb)
                print("[ERROR] Something is wrong! Please check the input tsv file.")
                exit(1)
            else:
                # Prepare array with cell type names
                cell_type_names = [ Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, tB ]

                # Prepare streams for output BED files and proceed with the next line of the tsv file
                cell_type_streams = []
                for cell_type_name in cell_type_names:
                    cell_type_streams.append(open(out_prefix + "_chicago_threshold_" + str(chicago_threshold) + "_" + cell_type_name + ".bed", 'wt'))
                continue

        # Check for trans-interactions and whether the baited digest comes before the other end in sequential order
        if dist == 'NA':
            # This is a trans-chromosomal interaction
            continue
        elif(int(float(dist)) < 0):
            # Note: The distances need to be converted to float first because some distances are given in scientific noation, e.g. '2e+05'
            sta = int(float(oeStart))
            end = int(float(baitEnd))
        else:
            sta = int(float(baitStart))
            end = int(float(oeEnd))

        interactions_coordinates = "chr" + baitChr + '\t' + str(sta) + '\t' + str(end)

        # Sanity check
        if sta > end:
            print(interactions_coordinates)
            print("[ERROR] Start coordiante is greater than end coordinate!")
            exit(1)

        # Get CHiCAGO scores for each cell type
        cell_type_scores = [Mon, Mac0, Mac1, Mac2, Neu, MK, EP, Ery, FoeT, nCD4, tCD4, aCD4, naCD4, nCD8, tCD8, nB, tB]
        for i in range(0,17):
            stream = cell_type_streams[i]
            score = float(cell_type_scores[i])
            if chicago_threshold <= score:
                #stream.write(interactions_coordinates + '\t' + str(score) + '\t' + cell_type_names[i] + '\n')
                stream.write(interactions_coordinates + '\n')

fp.close()

### Close streams for output files
##################################

for stream in cell_type_streams:
    stream.close()

print("[INFO] ... done.")