from collections import defaultdict

class ParseCallback:
    """
    Callback function superclass. This is intended to be used to extract
    specific data from the EnhancedInteractionFileParser
    """
    def __init__(self):
        super().__init__()
    
    def on_end_line(self,
                    chrompos_a, 
                    chrompos_b, 
                    syms_a, 
                    tsss_a, 
                    syms_b, 
                    tsss_b, 
                    enrichment_pair_tag, 
                    strand_pair_tag, 
                    interaction_category, 
                    neg_log_p_value, 
                    rp_simple, 
                    rp_twisted, 
                    i_dist):
        """
        function signature to parse elements from one line of the EnhancedInteraction file.
        """
        pass # todo consider data Q/C here

class Digest2SimpleTwistedCallback:
    """
    Create a dictionary of ChromosomalPosition to counts for all digests in the EnhancedInteraction file.
    """
    def __init__(self):
        super().__init__()
        self.stdict = defaultdict(dict)

    def on_end_line(self,
                    chrompos_a, 
                    chrompos_b, 
                    syms_a, 
                    tsss_a, 
                    syms_b, 
                    tsss_b, 
                    enrichment_pair_tag, 
                    strand_pair_tag, 
                    interaction_category, 
                    neg_log_p_value, 
                    rp_simple, 
                    rp_twisted, 
                    i_dist):
        self.stdict[chrompos_a] = {"S":rp_simple, "T": rp_twisted}
        self.stdict[chrompos_b] = {"S":rp_simple, "T": rp_twisted}

    def get_simple_twisted_dict(self):
        return self.stdict


class ShellDump:
    """
    Class to dump the parse result to the shell
    """
    def __init__(self):
        super().__init__()

    
    def on_end_line(self,
                    chrompos_a, 
                    chrompos_b, 
                    syms_a, 
                    tsss_a, 
                    syms_b, 
                    tsss_b, 
                    enrichment_pair_tag, 
                    strand_pair_tag, 
                    interaction_category, 
                    neg_log_p_value, 
                    rp_simple, 
                    rp_twisted, 
                    i_dist):
        print("chrompos_a {}, syms_a {} tsss_a {}".format(chrompos_a, syms_a, tsss_a))
        print("chrompos_b {}, syms_b {} tsss_b {}".format(chrompos_b, syms_b, tsss_b))
        print("enrichment_pair_tag {} strand_pair_tag {} interaction_category {}".format(enrichment_pair_tag, strand_pair_tag, interaction_category))
        print("neg_log_p_value {}, rp_simple {}, rp_twisted {}, i_dist {}".format(neg_log_p_value, rp_simple, rp_twisted, i_dist ))