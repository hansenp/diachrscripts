#!/usr/bin/env python

"""
This script takes a file in Diachromatic interaction format and a matching digest file exported from GOPHER and
generates a BedGraph file with an interaction profile.
"""

import argparse
import sys
from diachr import DiachromaticInteractionSet

# Parse command line
####################

parser = argparse.ArgumentParser(description='Create a BedGraph file with an interaction profile. The bins of the'
                                             'BEdGraph correspond to the restriction fragments.')
parser.add_argument('-o', '--out-prefix',
                    help='Common prefix for all generated files, which can also contain the path.',
                    default='OUT_PREFIX')
parser.add_argument('-i', '--diachromatic-interaction-file',
                    help='File in Diachromatic interaction format.',
                    required=True)
parser.add_argument('-s', '--select-top-n',
                    help='Use \'<N>:RPNUM\' or \'<N>:RPMAX\' to select the top n interactions sorted with respect to '
                         'total or maximum read pair count.',
                    default='ALL')
parser.add_argument('-ic', '--interaction-category',
                    help='Tag for interaction category.',
                    choices=['DIX', 'DI', 'UIR', 'UI', 'ALL'],
                    default='ALL')
parser.add_argument('-ec', '--enrichment-category',
                    help='Tag for enrichment category of interactions.',
                    choices=['NN', 'NE', 'EN', 'EE'],
                    default='ALL')
parser.add_argument('-t', '--profile-type',
                    help='Type of profile that will be created.'
                         'By default, an interaction profile (\'INUM\') is created.'
                         'With this type, the count in a bin for a given restriction fragment corresponds to the number'
                         'of interactions spanning that fragment.'
                         'For type \'RPNUM\', the bin count of a fragment corresponds to the total number of spanning'
                         'supporting read pairs, and for type \'RPMAX\','
                         'not the total number, but only the maximum of the counts for the four read types is taken'
                         'into account for each interaction.',
                    choices=['INUM', 'RPNUM', 'RPMAX'],
                    default='INUM')
parser.add_argument('-d', '--gopher-digest-file',
                    help='Digest file created with GOPHER using the same genome build and restriction enzyme that'
                         'were used to call the interactions.',
                    required=True)
parser.add_argument('-n', '--track-name',
                    help='Will be used as BedGraph track name and description.',
                    default='TRACK_NAME')
parser.add_argument('-c', '--track-color',
                    help='Will be used as BedGraph track color.',
                    default='0,0,0')
parser.add_argument('-b', '--write-bed', help='If true, an additional BED file will be created or development purposes'
                                              'containing the spanned regions of all interactions.'
                                              'Use this option only with small interaction datasets!'
                                              'If there are many interactions (100,000<), the BED file quickly becomes'
                                              'very large.',
                    default=False, action='store_true')
parser.add_argument('-v', '--verbose', help='If true, the script will provide feedback on progress.',
                    default=False, action='store_true')

# Assign command line arguments to script variables
args = parser.parse_args()
OUT_PREFIX = args.out_prefix
INTERACTION_FILE = args.diachromatic_interaction_file
SELECT_TOP_N = args.select_top_n
I_CAT = args.interaction_category
E_CAT = args.enrichment_category
PROFILE_TYPE = args.profile_type
DIGEST_FILE = args.gopher_digest_file
TRACK_NAME = args.track_name
TRACK_COLOR = args.track_color
WRITE_BED = args.write_bed
VERBOSE = args.verbose

# Define names for output files
BEDGRAPH_FILENAME = OUT_PREFIX + '_profile.bedgraph'
BED_FILENAME = OUT_PREFIX + '_interactions.bed'

# Report on parameters
if VERBOSE:
    parameter_info = "[INFO] " + "Input parameters" + '\n'
    parameter_info += "\t[INFO] --out-prefix: " + OUT_PREFIX + '\n'
    parameter_info += "\t[INFO] --diachromatic-interaction-file:" + '\n'
    parameter_info += "\t\t[INFO] " + INTERACTION_FILE + '\n'
    parameter_info += "\t[INFO] --diachromatic-interaction-file:" + '\n'
    parameter_info += "\t[INFO] --select-top-n: " + SELECT_TOP_N + '\n'
    parameter_info += "\t[INFO] --interaction-category: " + I_CAT + '\n'
    parameter_info += "\t[INFO] --enrichment-category: " + E_CAT + '\n'
    parameter_info += "\t[INFO] --profile-type: " + PROFILE_TYPE + '\n'
    parameter_info += "\t[INFO] --gopher-digest-file: " + DIGEST_FILE + '\n'
    parameter_info += "\t[INFO] --track-name: " + TRACK_NAME + '\n'
    parameter_info += "\t[INFO] --track-color: " + TRACK_COLOR + '\n'
    parameter_info += "\t[INFO] --write-bed: " + str(WRITE_BED) + '\n'
    parameter_info += "\t[INFO] --verbose: " + str(VERBOSE) + '\n'
    print(parameter_info + '\n')

# Prepare data structures
#########################

# Load interaction set
d11_interaction_set = DiachromaticInteractionSet(rpc_rule='ht')
d11_interaction_set.parse_file(
    i_file=INTERACTION_FILE,
    verbose=VERBOSE)

# Select the top n interactions sorted with respect to 'RPNUM' or 'RPMAX'
if SELECT_TOP_N != 'ALL':
    if VERBOSE:
        print()
        print('[INFO] Selecting top n interactions...')
    n, sort_criterion = SELECT_TOP_N.split(':')
    if sort_criterion not in ['RPNUM', 'RPMAX']:
        print('[ERROR] Invalid sort criterion! Must be \'RPNUM\' or \'RPMAX\'.')
    print('\t[INFO] Sorting parameters: ' + SELECT_TOP_N)
    d11_interaction_set_top_selected = \
        d11_interaction_set.sort_and_select_top_n_interactions(sort_by=sort_criterion, top_n=int(n))
    d11_interaction_set = d11_interaction_set_top_selected
    if VERBOSE:
        print('[INFO] ...done.')

if VERBOSE:
    print()
    print('[INFO] Reading digest file into data structure...')

D_MAP = {}
with open(DIGEST_FILE, 'rt') as fh:
    next(fh)
    for line in fh:
        fields = line.rstrip().split('\t')
        key = fields[0] + ':' + fields[1]
        D_MAP[key] = {'d_end': fields[2], 'cnt': 0}

if VERBOSE:
    print('[INFO] ...done.')

if VERBOSE:
    print()
    print('[INFO] Counting spanning interactions for all restriction fragments...')

# Write BED file for testing and development purposes
if WRITE_BED:
    fh_out_bed = open(BED_FILENAME, 'wt')
    track_name = TRACK_NAME + ' - BED test file'
    track_description = track_name
    fh_out_bed.write(
        'track type=bed name=\"' + track_name + '\" description=\"' + TRACK_NAME + '\" color=' + TRACK_COLOR + ' visibility=full' + '\n')

# Start counting
i_progress = 0
for d11_inter in d11_interaction_set.interaction_list:

    # Skip interactions
    if I_CAT != 'ALL' and I_CAT != d11_inter.get_category():
        continue
    if E_CAT != 'ALL' and E_CAT != d11_inter.enrichment_status_tag_pair:
        continue
    if d11_inter.chrA != d11_inter.chrB:
        continue

    # Define increment
    if PROFILE_TYPE == 'INUM':
        increment_by = 1
    elif PROFILE_TYPE == 'RPNUM':
        increment_by = d11_inter.rp_total
    elif PROFILE_TYPE == 'RPMAX':
        increment_by = max(d11_inter._simple_1, d11_inter._simple_2, d11_inter._twisted_1, d11_inter._twisted_2)
    else:
        print('[ERROR] Invalid increment type!')
        exit(0)

    # Increment all restriction fragments spanned by this interaction
    chrom = d11_inter.chrA
    i_sta = d11_inter.fromA
    i_end = d11_inter.toB
    d_sta = i_sta
    key = chrom + ':' + str(d_sta)
    while key in D_MAP and int(d_sta) < i_end:
        D_MAP[key]['cnt'] += increment_by
        d_sta = int(D_MAP[key]['d_end']) + 1
        key = chrom + ':' + str(d_sta)

    # Write interaction to BED file for testing and development purposes
    if WRITE_BED:
        name = d11_inter.get_category() + '|' + d11_inter.enrichment_status_tag_pair
        fh_out_bed.write(chrom + '\t' + str(i_sta - 1) + '\t' + str(i_end) + '\t' + name + '\n')

    # Provide feedback on progress
    if VERBOSE:
        if i_progress > 0 == i_progress % 100000:
            print("\t[INFO] Parsed " + "{:,}".format(i_progress) + " interaction lines ...")
        i_progress += 1

if WRITE_BED:
    fh_out_bed.close()

if VERBOSE:
    print('[INFO] ...done.')

# Write BedGraph file with profile
##################################

if VERBOSE:
    print()
    print('[INFO] Write BedGraph file with profile...')

bedgraph_fh = open(BEDGRAPH_FILENAME, 'wt')
bedgraph_fh.write(
    'track type=bedGraph name=\"' + TRACK_NAME + '\" description=\"' + TRACK_NAME + '\" color=' + TRACK_COLOR + ' autoScale=on alwaysZero=on visibility=full priority=20' + '\n')
with open(DIGEST_FILE, 'rt') as fp:
    next(fp)
    for line in fp:
        fields = line.rstrip().split('\t')
        key = fields[0] + ':' + fields[1]
        chrom = key.split(':')[0]
        sta = key.split(':')[1]
        end = D_MAP[key]['d_end']
        cnt = D_MAP[key]['cnt']
        bedgraph_fh.write(chrom + '\t' + str(int(sta) - 1) + '\t' + end + '\t' + str(cnt) + '\n')
bedgraph_fh.close()

if VERBOSE:
    print('[INFO] ...done.')
    print()
    generated_file_info = "[INFO] Generated files:" + '\n'
    generated_file_info += "\t[INFO] " + BEDGRAPH_FILENAME + '\n'
    if WRITE_BED:
        generated_file_info += "\t[INFO] " + BED_FILENAME + '\n'
    print(generated_file_info + '\n')
