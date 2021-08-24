from diachr.diachromatic_interaction import DiachromaticInteraction
from typing import Tuple
import os
from collections import defaultdict


class TadRegion:
    """
    simple class with start and end positions of TADs
    
    """
    def __init__(self, chr:str, s:int, e: int) -> None:
        self._chr = chr
        self._start = s
        self._end = e

    def contains(self, pos: int):
        """
        return true if this TAD contains the indicated position
        Note that we assume we are on the same chromosome
        """
        return self._start <= pos <= self._end

    

class TadRegion:

    def __init__(self, tad_input: str):
        if not os.path.isfile(tad_input):
            raise FileNotFoundError("Could not find TAD BED file {}".format(tad_input))
        chr_d = defaultdict(list)
        with open(tad_input) as f:
            next(f) # discard header
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 3:
                    raise ValueError("Malformed BED3 line: {}".format(line))
                chr = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                tadpos = TadRegion(chr=chr, s=start, e=end)
                chr_d[chr].append(tadpos)  
        print("[INFO] Parsed {} chromosomes from {}".format(len(chr_d), tad_input))   
        # Chromosomes should be sorted from the input file, but we will
        # sort here to be sure
        for chrom, tadlist in chr_d.items():
            # ut.sort(key=lambda x: x.count, reverse=True)
            tadlist.sort(key=lambda x: x._start)
        self._chr_d = chr_d
        
   



    def get_distance_to_neighbouring_boundaries(chr:str, start:int, end:int) -> Tuple[int,int]:
        """
        Get distance to the closest TAD boundaries
        If the argument is located within a boundary, what happens?
        """
        pass


    def get_overlapped_boundary_count_from_interaction(self, ia: DiachromaticInteraction) -> int:
        chrA = ia.chrA
        chrB = ia.chrB
        if chrA != chrB:
            raise ValueError("Error -- two different chromosomes")
        start = ia.fromA
        end = ia.toB
        return self.get_overlapped_boundary_count(chr=chrA, start=start, end=end)
        


    def get_overlapped_boundary_count(self, chr:str, start:int, end:int) -> int:
        """
        The argument is typically an interaction. start is the 'left' digest and end is the 'right' digest
        """
        if not chr in self._chr_d:
            raise ValueError("Could not find chr ({}) in chrom dict".format(chr))
        tadlist = self._chr_d[chr]
        i = -1
        j = -1
        m = 0
        for tad in tadlist:
            if tad.contains(start):
                i = m
            if tad.contains(end):
                j = m
            m += 1
        if i < 0 or j < 0:
            raise ValueError("COULD NOT FIND IT")
        return 1 + j - i
