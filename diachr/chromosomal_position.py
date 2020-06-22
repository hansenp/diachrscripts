


class ChromosomalPosition:
    """
    Class to coordinate use of chromosomal positions such as chr1:123-456
    """
    def __init__(self, coordinate_string):
        if not ":" in coordinate_string:
            raise ValueError("Chromosomal coordinate string did not contain semicolon")
        if not "-" in coordinate_string:
            raise ValueError("Chromosomal coordinate string did not contain dash")
         # Extract digest coordinates
        fields = coordinate_string.split(":")
        if len(fields) != 2:
            raise ValueError("Chromosomal coordinate string had bad split on :")
        self._chromosome = fields[0]
        pos = fields[1].split("-")
        if len(pos) != 2:
            raise ValueError("Chromosomal coordinate string had bad split on -")
        self._start = int(pos[0])
        self._end = int(pos[1])

    @property
    def chromosome(self):
        return self._chromosome
    
    @property
    def start(self):
        return self._start
    
    @property
    def end(self):
        return self._end

    def __hash__(self):
        return hash(self.chromosome) + 31*hash(self.start) + 31*hash(self.end)

    def __eq__(self, other):
        return isinstance(other, self.__class__) and other.chromosome == self.chromosome and other.start == self.start and other.end == self.end

    def __str__(self):
        return "%s:%d-%d" %  (self.chromosome, self.start, self.end)