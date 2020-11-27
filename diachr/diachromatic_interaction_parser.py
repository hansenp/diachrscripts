import gzip
import os
from typing import Tuple, List
from collections import defaultdict
from .diachromatic_interaction import DiachromaticInteraction

class DiachromaticInteractionParser:
    """
    This class coordinates the parsing of Diachromatic interaction files.
    It creates a list of EnhancedInteraction objects, each representing one data line.
    """

    def __init__(self):
        self._file_num = 0
        self._inter_dict = defaultdict(DiachromaticInteraction)

    def _parse_line(self, line: str) -> DiachromaticInteraction:
        F = line.rstrip().split('\t')
        if len(F) < 9:
            raise ValueError("Malformed line {}".format(line))

        chrA = F[0]
        fromA = int(F[1])
        toA = int(F[2])
        statusA = F[3]
        chrB = F[4]
        fromB = int(F[5])
        toB = int(F[6])
        statusB = F[7]
        st_string = F[8]  # something like 2:1, representing S:T, simple and twisted counts
        st_fields = st_string.split(":")
        if len(st_fields) != 2:
            raise ValueError("Malformed simple:twisted field in diachromatic line: " + line)
        simple = int(st_fields[0])
        twisted = int(st_fields[1])

        di_inter = DiachromaticInteraction(chrA=chrA, fromA=fromA, toA=toA, statusA=statusA,
                                          chrB=chrB, fromB=fromB, toB=toB, statusB=statusB, simple=simple,
                                          twisted=twisted)
        return di_inter

    def parse_file(self, path) -> List:

        if not os.path.exists(path):
            raise ValueError("Could not find file at %s" % path)

        self._file_num += 1

        if path.endswith(".gz"):
            with gzip.open(path, 'rt') as fp:
                n_progress = 0
                for line in fp:
                    n_progress += 1
                    if n_progress % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_progress) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.simple, twisted=d_inter.twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter
        else:
            with open(path) as fp:
                n_progress = 0
                for line in fp:
                    n_progress += 1
                    if n_progress % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_progress) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.simple, twisted=d_inter.twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter

    @property
    def i_list(self):
        return self._inter_dict.values()

    @property
    def file_num(self):
        return self._file_num
