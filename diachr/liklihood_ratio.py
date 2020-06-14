from scipy.stats.distributions import chi2
from scipy.stats import binom
from math import log


class DirectionalityLikelihoodRatio:
    """
    A class to calculate the likelihood ratio for there being a different
    directionality at a particular digest. The input consists of two
    directionality dictionaries as produced from two extended interaction
    files from the EnhancedInteractionFileParser using the callback
    Digest2SimpleTwistedCallback
    """
    def __init__(self, experiment_a, experiment_b):
        if not isinstance(experiment_a, dict):
            raise TypeError("experiment_a must be a dictionary!")
        if not isinstance(experiment_b, dict):
            raise TypeError("experiment_b must be a dictionary!")
        self.experiment_a = experiment_a
        self.experiment_b = experiment_b

    def _get_combined_log_likelihood(self, simple_a, twisted_a, simple_b, twisted_b):
        total = simple_a + twisted_a + simple_b + twisted_b
        total_simple = simple_a + simple_b
        p = total_simple/total
        k = total_simple
        n = total
        return binom.logpmf(k, n, p)

    def _get_separated_log_likelihood(self, simple_a, twisted_a, simple_b, twisted_b):
        total_a = simple_a + twisted_a
        total_b = simple_b + twisted_b
        p_a = simple_a / total_a
        p_b = simple_b / total_b
        k = simple_a
        n = total_a
        pmf_a =  binom.logpmf(k, n, p_a)
        k = simple_b
        n = total_b
        pmf_b =  binom.logpmf(k, n, p_b)
        return pmf_a + pmf_b




    def compare_interactions(self):
        print("[INFO] Calculating likelihood ratios for A (n=%d) and B (n=%d)" %
            (len(self.experiment_a), len(self.experiment_b)))
        for k, v in self.experiment_a.items():
            if not k in self.experiment_b:
                continue
            simple_a = v.get("S")
            twisted_a = v.get("T")
           
            b = self.experiment_b[k]
            simple_b = b.get("S")
            twisted_b = b.get("T")
           
            ll_combined = self._get_combined_log_likelihood(simple_a, twisted_a, simple_b, twisted_b)
            ll_separated = self._get_separated_log_likelihood(simple_a, twisted_a, simple_b, twisted_b)
            LR = 2*(ll_combined-ll_separated)
            p = chi2.sf(LR, 1)
            if p < 0.02:
                print(k, p)

