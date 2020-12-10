from unittest import TestCase
import os
import sys

from numpy.lib.function_base import place
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from diachr import BinomialModel

class TestBinomial(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.binomial_model = BinomialModel()

    def test_binomial_pval_18_7(self):
        """
        In R
        binom.test(18, 25, 0.5, alternative="two.sided") 
        number of successes = 18, number of trials = 25, p-value = 0.04329
        This is effectively the same as
        2*(1-pbinom(17,25,0.5)) 
        and calculates the probability of having scores ranging from 0 to 7 and 18 to 25 (two-sided test). 
        """
        simple = 18
        twisted = 7
        pval = self.binomial_model.get_binomial_p_value(n_simple=simple, n_twisted=twisted)
        self.assertAlmostEqual(0.04329, pval, places=5) # equal up to 5 significant digits


    def test_binomial_pval_2_34(self):
        """
        In R
        binom.test(2, 36, 0.5, alternative="two.sided") 
        number of successes = 2, number of trials = 36, p-value = 1.941e-08
        """
        simple = 2
        twisted = 34
        pval = self.binomial_model.get_binomial_p_value(n_simple=simple, n_twisted=twisted)
        self.assertAlmostEqual(1.941e-08, pval, places=5) # equal up to 5 significant digits