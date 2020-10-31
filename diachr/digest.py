
#######################################################################################################################

class Digest:

    # Class to represent a genomic region that corresponds to a restriction digest

    # Attributes
    chromosome = None
    start = None
    end = None
    active = False

    # Initializer

    def __init__(self, chromosome:str, start:int, end:int, active:bool = False):
        self._chromosome = chromosome
        self._start = start
        self._end = end
        self._active = active

    # Methods
    def set_active(self):
        self._active = True

    def is_active(self):
        return self._active

    def get_chromosome(self):
        return self._chromosome

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    @property
    def coordinate_string(self):
        return self._chromosome + ":" + str(self._start) + "-" + str(self._end)

##############################