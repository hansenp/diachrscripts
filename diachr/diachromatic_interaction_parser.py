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
        self._file_dict = {}
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

    def parse_file(self, i_file):
        """
        Parses a file with interactions. For interactions that have already been parsed from another file,
        only the simple and twisted read pair counts are added.
        :param i_file: File with interactions
        """
        if not os.path.exists(i_file):
            raise ValueError("Could not find file at %s" % i_file)

        if i_file.endswith(".gz"):
            with gzip.open(i_file, 'rt') as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if n_lines % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.simple, twisted=d_inter.twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter
        else:
            with open(i_file) as fp:
                n_lines = 0
                for line in fp:
                    n_lines += 1
                    if n_lines % 1000000 == 0:
                        print("\t[INFO] Parsed " + str(n_lines) + " interaction lines ...")
                    d_inter = self._parse_line(line)
                    if d_inter.key in self._inter_dict:
                        self._inter_dict[d_inter.key].append_interaction_data(simple=d_inter.simple, twisted=d_inter.twisted)
                    else:
                        self._inter_dict[d_inter.key] = d_inter

        self._file_dict[i_file] = n_lines

    def get_file_dict_info(self):
        """
        :return: String that contains information about parsed data.
        """
        file_dict_info = "[INFO] Read interaction data from " + str(len(self._file_dict)) + " files:" + '\n'
        for k, v in self._file_dict.items():
            file_dict_info += "\t[INFO] " + str(v) + " interactions from " + k  + '\n'
        file_dict_info += "[INFO] The union of all interactions has "  + str(len(self._inter_dict)) + " interactions." + '\n'
        return(file_dict_info)

    def write_diachromatic_interaction_file(self, required_replicates = 1, target_file_name = None):
        """
        Writes interactions that occur in a specified minimum number of replicates to a file
        and returns a string with information on this writing process.

        :param required_replicates: Minimum number of replicates
        :param target_file_name: Generated file with interactions
        :return: String with information on this writing process
        """
        write_info = "[INFO] Writing interactions that occur in at least " + str(required_replicates) + " replicates to " + target_file_name + '\n'
        out_fh = gzip.open(target_file_name, 'wt')
        n_has_required_data = 0
        n_incomplete_data = 0
        for d_inter in self._inter_dict.values():
            if d_inter.has_data_for_required_replicate_num(required_replicates):
                n_has_required_data += 1
                out_fh.write(d_inter.get_diachromatic_interaction_line() + '\n')
            else:
                n_incomplete_data += 1
        out_fh.close()
        write_info += "\t[INFO] Interactions that occur in at least " + str(required_replicates) + " replicates: " + str(n_has_required_data) + '\n'
        write_info += "\t[INFO] Other interactions: " + str(n_incomplete_data) + '\n'
        write_info += '\n'
        write_info += "FILE_NAME" + "\t" + "INTERACTIONS_NUMBERS" + "\t" + "REQUIRED_INTERACTIONS" + "\t" + "HAS_ALL_DATA" + "\t" + "INCOMPLETE_DATA" + '\n'
        write_info += str(target_file_name) + "\t" + str(list(self._file_dict.values())) + "\t" + str(required_replicates) + "\t" + str(n_has_required_data) + "\t" + str(n_incomplete_data) + '\n'
        return write_info

    @property
    def file_num(self):
        return len(self._file_dict)
