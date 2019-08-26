import diachrscripts_toolkit as dclass
from scipy.stats import binom
import pytest
from math import log, log10

#pytestmark = pytest.mark.filterwarnings("ignore")

class TestKinteractionPvalue:

    def test_compare_scipi_sf_and_one_minus_cdf_pvalue(self):
        """
        Check whether survival function corresponds to P-value.

        For small k, n and x, SF(k,n,p=2*(0.5**x)) == 1-CDF(k,n,p=2*(0.5**x)) but for larger values the term with
        CDF evaluates to zero (presumably when CDF(k,n,p=2*(0.5**x) becomes smaller than the machine epsilon).
        Therefore, we use 'logsf' in order to calculate the P-value.

        Note: This test will fail for larger k, n and x.
        """
        k = 2
        n = 10
        x = 3
        p = 2 * (0.5 ** x)
        my_logsf_p_val, error_state = dclass.get_x_inter_binomial_log10_p_value(k, n, x)
        assert my_logsf_p_val == -log10(1-binom.cdf(k,n,p))

