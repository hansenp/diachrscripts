import os
import gzip
from collections import defaultdict
from .diachromatic_interaction import DiachromaticInteraction



class DiachromaticParser:
    """
    Parses Diachromatic interaction files (https://github.com/TheJacksonLaboratory/diachromatic).
    We assume that multiple experiments have been done with the same biological samples and that
    this script will combine these files into a single output file intended for downstream analysis.
    The constructor takes the path to a directory where these files are located and assumes that
    all files with the suffix .tsv.gz are Diachromatic interaction files.
    """
    def __init__(self, path:str, required_replicates: int = None, prefix: str = "INTERACTIONS") -> None:
        """
        Parses one or more Diachromatic interaction files
        Args:
        path (str): Path to a directory with interaction files or to a single interaction file.
        required_replicates (int): Minimum number of interactions to be included in a combined interaction file. Required if a directory is passed, otherwise not relevant
        prefix (str): prefix for outfiles
        """
        self._path = path
        self._gzfiles = []
        self._required_replicates = required_replicates
        self._prefix = prefix
        self._interaction_d = defaultdict(DiachromaticInteraction)
        self._n_iteraction = [] # list with total number of interactions per file
        if os.path.isdir(self._path):
            print("[INFO] Parsing diachromatic files in  directory at ", self._dir)
            if not isinstance(required_replicates, int):
                raise ValueError("Need to pass an integer value for required_replicates")
            for file in os.listdir(self._path):
                if file.endswith(".tsv.gz"):
                    gzpath = os.path.join(self._path, file)
                    self._gzfiles.append(gzpath)
            if len(self._gzfiles) < self._required_replicates:
                raise ValueError("Not enough replicates. We found only %d diachromatic files be we require %d" % (len(self._gzfiles), self._required_replicates))
            print("[INFO] We found %d Diachromatic interaction files at %s" % (len(self._gzfiles), self._path))         
            for gzfile in self._gzfiles:
                self._parse_gzip_file(gzfile=gzfile)
        elif os.path.isfile(self._path):
            self._parse_gzip_file(self._path)
        else:
            raise ValueError("path ({}) must point to a directory or to a file".format(path))
            

    def _parse_gzip_file(self, gzfile: str):
        """
        Parse a single gzip Diachromatic file.
        Note that this can be either the original file from Diachromatic from 
        a single experiment, or it can be the combined file
        # format -
        # chr5\t169261149\t169264959\tI\tchr5\t169981022\t169985681\tI\t0:2
        """
        if not os.path.isfile(gzfile):
            raise ValueError("We were expected a single gzip diachromatic file")
        n_iteraction = 0
        with gzip.open(gzfile, 'rt') as f:
            for line in f:
                n_iteraction += 1
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
                st_string = F[8] # something like 2:1, representing S:T, simple and twisted counts
                st_fields = st_string.split(":")
                if len(st_fields) != 2:
                    raise ValueError("Malformed simple:twisted field in diachromatic line: " + line)
                simple = int(st_fields[0])
                twisted = int(st_fields[1])
                iaction = DiachromaticInteraction(chrA=chrA, fromA=fromA, toA=toA, statusA=statusA,
                        chrB=chrB, fromB=fromB, toB=toB, statusB=statusB, simple=simple, twisted=twisted)
                if iaction.key in self._interaction_d:
                    self._interaction_d[iaction.key].append_interaction_data(simple=simple, twisted=twisted)
                else:
                    self._interaction_d[iaction.key] = iaction
        self._n_iteraction.append(n_iteraction)

    def get_interaction_dict(self):
        return self._interaction_d


    def write_combined_interaction_file(self):
        fname = self._prefix + "_at_least_in_" + str(self._required_replicates) + "_replicates_interactions.tsv.gz"
        outfh = gzip.open(fname, 'wt')
        print("[INFO] Writing interactions to file ...")
        n_has_required_data = 0
        n_incomplete_data = 0
        for _, iaction in self._interaction_d.items():
            if not iaction.has_data_for_required_replicate_num(self._required_replicates):
                n_incomplete_data += 1
            else:
                n_has_required_data += 1
                outfh.write(iaction.output_summary() + "\n")

            if (n_incomplete_data + n_has_required_data)%1000000==0:
                print("\t[INFO] " + (str(n_incomplete_data + n_has_required_data)) + " interactions processed ...")
        outfh.close()
        print("[INFO] We wrote all interactions to file: {}".format(fname))
        fname = self._prefix + "_at_least_in_" + str(self._required_replicates) + "_replicates_summary.txt"
        outfh = open(fname, 'wt')
        outfh.write("OUT_PREFIX" + "\t" + "INTERACTIONS_NUMBERS" + "\t" + "REQUIRED_INTERACTIONS" + "\t" + "HAS_ALL_DATA" + "\t" + "INCOMPLETE_DATA" + "\n")
        outfh.write(self._prefix  + "\t" + str(self._n_iteraction) + "\t" + str(self._required_replicates) + "\t" + str(n_has_required_data) + "\t" + str(n_incomplete_data) + "\n")
        outfh.close()

        print("[INFO] Interactions with at least {} data points: {}, lacking interactions: {}".format(self._required_replicates, n_has_required_data, n_incomplete_data))
        print("[INFO] We wrote summary statistics to file: {}".format(fname))
        print("[INFO] ... done.")

    

    
