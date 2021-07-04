import argparse
import os
from diachr import DirectedTadVisualizer

parser = argparse.ArgumentParser(description='Create directed TAD visualization.')
parser.add_argument('--infile','-i', required=True,  help='path to directed interaction file')
parser.add_argument('--chrom','-c',required=True,help='chromosome')
parser.add_argument('--begin','-b',required=True,help='begin position')
parser.add_argument('--end','-e',required=True,help='end position')
args = parser.parse_args()
infile = args.infile
chrom = args.chrom
begin = int(args.begin)
end = int(args.end)
print("[INFO] Creating graphic for {}".format(infile))
print("[INFO] Chromosome {}:{}-{}".format(chrom, begin, end))

visualizer = DirectedTadVisualizer(fname=infile)
outname = 'dtadvis-out.tsv'
if not os.path.isfile(outname):
    visualizer.extract_interactions(chrom=chrom, begin=begin, end=end, fname=outname)
visualizer.create_visualization(outname)

