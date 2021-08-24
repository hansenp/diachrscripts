from diachr.diachromatic_interaction import DiachromaticInteraction
from typing import Tuple
import os
from collections import defaultdict



class TadBoundary:
    def __init__(self, chr:str, b:int, e:int)->None:
        self._chr = chr
        self._start = b
        self._end = e



class TadBoundarySet:
    def __init__(self, input_bedfile: str)->None:
        if not os.path.isfile(input_bedfile):
            raise FileNotFoundError("Could not find TAD BED file {}".format(input_bedfile))
        chr_d = defaultdict(list)
        with open(input_bedfile) as f:
            next(f) # discard header
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 3:
                    raise ValueError("Malformed BED3 line: {}".format(line))
                chr = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                tboundary = TadBoundary(chr=chr, b=start, e=end)
                chr_d[chr].append(tboundary)