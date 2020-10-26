from typing import Tuple, List


class EnhancedInteraction:
    """
    chr_a: Chromosome of the first digest of the interaction (e.g. chr1)
    sta_a: Start coordinate of the first digest of the interaction (e.g. 123456)
    end_a: End coordinate of the first digest of the interaction (e.g. 234567)
    syms_a: Comma separated list of gene symbols associated with TSS on the first digest (e.g. HIST1H1E,HIST1H2BD)
    tsss_a: Comma separated list of TSS coordinates on the first digest (e.g. chr6:26158121:+,chr6:26156331:+)
    chr_b, sta_b, end_b, syms_b, tsss_b: Same as for the first digest but for the second digest
    enrichment_pair_tag: Two letter tag indicating the enrichment status of the two digests (e.g.  AA, AI, II)
    strand_pair_tag: Two symbol tag separated by '/' indicating the strands of TSS on first and second digest (e.g. -/-, +/+, -/+, +/-, -/d, ...)
    interaction_category: Interaction category wi
    """
    def __init__(self, chrA:str, staA:int, endA:int, symsA:str, tsssA:str,chrB:str, staB:int, endB:int, symsB:str, tsssB:str,
                 enr_pair_tag:str, strand_pair_tag:str, intera_cat:str, neg_log_p:float, read_pair_left:int, read_pair_right:int,i_dist:int):
        self._chr_a = chrA
        self._start_a = staA
        self._end_a = endA
        self._syms_a = symsA
        self._symbols_a = symsA.split(",")
        self._tsss_a = tsssA
        self._trans_starts_a = tsssA.split(",")
        self._chr_b = chrB
        self._start_b = staB
        self._end_b = endB
        self._syms_b = symsB
        self._symbols_b = symsB.split(",")
        self._tsss_b = tsssB
        self._trans_starts_b = tsssB.split(",")
        self._enrichment_pair_tag = enr_pair_tag
        self._strand_pair_tag = strand_pair_tag
        self._interaction_category = intera_cat
        self._neg_log_p_value = neg_log_p
        self._read_pair_left = read_pair_left
        self._read_pair_right = read_pair_right
        self._i_dist = i_dist

    @property
    def chr_a(self) -> str:
        return self._chr_a

    @property
    def sta_a(self) -> int:
        return self._start_a

    @property
    def end_a(self) -> int:
        return self._end_a

    @property
    def chr_b(self) -> str:
        return self._chr_b

    @property
    def sta_b(self) -> int:
        return self._start_b

    @property
    def end_b(self) -> int:
        return self._end_b

    @property
    def syms_a(self) -> str:
        return self._syms_a

    @property
    def syms_b(self) -> str:
        return self._syms_b

    @property
    def tsss_a(self) -> str:
        return self._tsss_a

    @property
    def tsss_b(self) -> str:
        return self._tsss_b

    @property
    def enrichment_pair_tag(self) -> str:
        return self._enrichment_pair_tag

    @property
    def strand_pair_tag(self):
        return self._strand_pair_tag

    @property
    def interaction_category(self) -> str:
        return self._interaction_category

    @property
    def neg_log_p_value(self) -> float:
        return self._neg_log_p_value

    @property
    def read_pair_left(self) -> int:
        return self._read_pair_left

    @property
    def read_pair_right(self) -> int:
        return self._read_pair_right

    @property
    def i_dist(self) -> int:
        return self._i_dist

    @property
    def rp_total(self) -> int:
        return self._read_pair_right + self._read_pair_left

    @property
    def coordinate_key_a(self):
        return "%s\t%d\t%d" % (self._chr_a, self._start_a, self._end_a)

    @property
    def coordinate_key_b(self):
        return "%s\t%d\t%d" % (self._chr_b, self._start_b, self._end_b)






