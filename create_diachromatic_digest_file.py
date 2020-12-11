#!/usr/bin/env python

"""
This script takes a BED file with coordinates of digests selected for enrichment as well as a Diachromatic digest file
and creates a new Diachromatic digest file in which digest regions that are contained in the BED file are marked with a
 '1' and all other digests are marked with '0'.

You can find a detailed documentation on this script in the relevant section in the RTD of this repository.
"""

import argparse


### Parse command line
######################

parser = argparse.ArgumentParser(description='Overwrite enrichment status in column 12 of a Diachromatic digest file.')
parser.add_argument('-o','--out-prefix', help='Prefix for output.', default='OUTPREFIX')
parser.add_argument('--enriched-digests-file', help='BED file with digests that were selected for target enrichment.', required=True)
parser.add_argument('-i','--diachromatic-digest-file', help='Diachromatic digest file created with GOPHER.', required=True)
parser.add_argument('-v','--verbose', help='Try to get more information about digests that were not found.', default=False, action='store_true')

args = parser.parse_args()
out_prefix = args.out_prefix
enriched_digests_file = args.enriched_digests_file
diachromatic_digest_file = args.diachromatic_digest_file
verbose = args.verbose

print("[INFO] " + "Input parameters")
print("\t[INFO] --out-prefix: " + out_prefix)
print("\t[INFO] --enriched-digests-file: " + str(enriched_digests_file))
print("\t[INFO] --diachromatic-digest-file: " + diachromatic_digest_file)


### Read BED file with coordinates of enriched digests
######################################################

print("[INFO] Reading list with digests selected for enrichment ...")

enriched_digests_set = set()
with open(enriched_digests_file, 'rt') as fp:
    for line in fp:
        chr, sta, end = line.rstrip().split('\t')
        enriched_digests_set.add(chr + '\t' + str(sta) + '\t' + str(end))

n_enriched_digests = len(enriched_digests_set)
print("\t[INFO] Read " + str(n_enriched_digests) + " digests ...")
print("[INFO] ... done.")


### Iterate input digest file and write new digest file with enrichment status flags overwritten
################################################################################################

print("[INFO] Iterating input digest file and writing new digest file ...")

# Create new Diachromatic digest file
out_diachromatic_digest_file_fn =  out_prefix + "_diachromatic_digest_file.txt"
out_diachromatic_digest_file_fh = open(out_diachromatic_digest_file_fn, 'wt')

enriched_digests_found_set = set()
enriched_digests_not_found_set = set()
with open(diachromatic_digest_file, 'rt') as fp:
    for line in fp:

        # Skip first line

        # Parse one line of a Diachromatic digest file
        fields = line.rstrip('\n').split('\t')
        chr = fields[0]        # Chromosome
        sta = fields[1]        # Fragment_Start_Position
        end = fields[2]        # Fragment_End_Position
        frN = fields[3]         # Fragment_Number
        rS5 = fields[4]         # 5'_Restriction_Site
        rS3 = fields[5]         # 3'_Restriction_Site
        Len = fields[6]         # Length
        gc5 = fields[7]         # 5'_GC_Content
        gc3 = fields[8]         # 3'_GC_Content
        rp5 = fields[9]         # 5'_Repeat_Content
        rp3 = fields[10]        # 3'_Repeat_Content
        selected = fields[11]  # Selected
        pr5 = fields[12]       # 5'_Probes
        pr3 = fields[13]       # 3'_Probes

        # Check whether digest is in enriched digest set
        SELECTED = 'F'
        PR5 = 0
        PR3 = 0
        key = chr + '\t' + sta + '\t' + end
        if key in enriched_digests_set:
            enriched_digests_found_set.add(key)
            SELECTED = 'T'
            PR5 = 1
            PR3 = 1

        # Write modified line to new digest file
        out_diachromatic_digest_file_fh.write(
            chr + '\t' +
            sta + '\t' +
            end + '\t' +
            frN + '\t' +
            rS5 + '\t' +
            rS3 + '\t' +
            Len + '\t' +
            gc5 + '\t' +
            gc3 + '\t' +
            rp5 + '\t' +
            rp3 + '\t' +
            SELECTED + '\t' +
            str(PR5) + '\t' +
            str(PR3) + '\n'
            )

out_diachromatic_digest_file_fh.close()

enriched_digests_not_found_set = enriched_digests_set - enriched_digests_found_set
n_enriched_digests_matched = len(enriched_digests_found_set)
n_enriched_digests_unmatched = len(enriched_digests_not_found_set)

print("\t[INFO] From the " + str(n_enriched_digests) + " enriched digests in the input file, "
      + str(n_enriched_digests_matched) + " were found in the digest file.")
print("\t[INFO] " + str(n_enriched_digests_unmatched) + " digest were not found.")

print("[INFO] ... done.")

# Write the coordinates of digests not found to a BED file
if 0 < n_enriched_digests_unmatched:
    out_digests_not_found_fn = out_prefix + "_digests_not_found.bed"
    out_digests_not_found_fh = open(out_digests_not_found_fn, 'wt')
    for d in enriched_digests_not_found_set:
        out_digests_not_found_fh.write(d + '\n')
    out_digests_not_found_fh.close()
    print("[INFO] Digests not found were written to: " + out_digests_not_found_fn)
    print("[INFO] Run the script again with the --verbose option to get more information about digests not found.")

if verbose == False:
    exit(0)


### Try to get more information about digests not found
#######################################################

print("[INFO] Trying to get more information about digests that were not found in the digest file ...")

add_range = 100
id_not_found = 0
for d in enriched_digests_not_found_set:
    id_not_found += 1
    print("------------------------------------------------------")
    print("DIGEST NOT FOUND " + str(id_not_found) + '\n')

    print("Digest coordinates:")
    print(d + '\n')

    CHR = d.split('\t')[0]
    STA = int(d.split('\t')[1])
    END = int(d.split('\t')[2])
    LEN = END - STA

    print("Digest length:")
    print(str(LEN) + '\n')

    with open(diachromatic_digest_file, 'rt') as fp:
        for line in fp:

            fields = line.rstrip('\n').split('\t')
            chr = fields[0]        # Chromosome
            sta = fields[1]        # Fragment_Start_Position
            end = fields[2]        # Fragment_End_Position

            if chr == CHR and (\
                    (STA - add_range <= int(sta) and int(sta) <= STA + add_range) or \
                    (END - add_range <= int(end) and int(end) <= END + add_range)):
                print("Associated line in input digest file:")
                print(line)

print("[INFO] ... done.")
