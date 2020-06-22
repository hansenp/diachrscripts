from diachr import EnhancedInteractionFileParser, Digest2SimpleTwistedCallback, ShellDump, DirectionalityLikelihoodRatio
# This script intends to ingest two extended analysis files and to compare digests for the percentage of directed interactions.




class Measurement:
    def __init__(self):
        super().__init__()

class Digest:
    def __init__(self):
        super().__init__()


acd4file = '/home/peter/data/di_and_uri_results/exact_ref/JAV_ACD4_RALT_0.0019_enhanced_interaction_file_with_di_and_uir.tsv.gz'
epfile = '/home/peter/data/di_and_uri_results/exact_ref/JAV_EP_RALT_0.0017_enhanced_interaction_file_with_di_and_uir.tsv.gz'

# Parse the ACD4 data
digest2simple = Digest2SimpleTwistedCallback()
parser = EnhancedInteractionFileParser(acd4file)
parser.parse([digest2simple])
acd4 = digest2simple.get_simple_twisted_dict()

# Parse the EP data
digest2simpleEp = Digest2SimpleTwistedCallback()
parser = EnhancedInteractionFileParser(epfile)
parser.parse([digest2simpleEp])
ep = digest2simpleEp.get_simple_twisted_dict()

likratio = DirectionalityLikelihoodRatio(acd4, ep)
likratio.compare_interactions()


