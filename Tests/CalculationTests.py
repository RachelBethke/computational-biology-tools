import unittest
import pandas as pd
import numpy as np
# Import your calculation functions
from Calculations import calc_pi, calc_watterson

class TestGenomicsCalculations(unittest.TestCase):
    
    def test_calc_pi_basic(self):
        """
        Test the calc_pi function with a basic haplotype block.
        """
        # Example haplotype block (4 individuals, 3 loci)
        data = {
            0: [0, 1, 0, 1],
            1: [1, 0, 1, 0],
            2: [0, 1, 1, 0]
        }
        haplotypes = pd.DataFrame(data)
        expected_pi = (0.333 + 0.333 + 0.167) / 3
        pi_result = calc_pi(haplotypes)
        self.assertAlmostEqual(pi_result, expected_pi, places=3)
    
    def test_calc_pi_empty(self):
        """
        Test calc_pi with an empty DataFrame (edge case).
        """
        empty_haplotypes = pd.DataFrame()
        self.assertEqual(calc_pi(empty_haplotypes), 0)
    
    def test_calc_pi_single_locus(self):
        """
        Test calc_pi with a DataFrame containing only one locus.
        """
        data = {0: [0, 1, 1, 0]}
        haplotypes = pd.DataFrame(data)
        expected_pi = 0.333
        pi_result = calc_pi(haplotypes)
        self.assertAlmostEqual(pi_result, expected_pi, places=3)
    
    # ---- TEST calc_watterson FUNCTION ----
    def test_calc_watterson(self):
        """
        Test the calc_watterson function with a basic haplotype block.
        """
        data = {
            0: [0, 1, 0, 1],
            1: [1, 0, 1, 0],
            2: [0, 1, 1, 0]
        }
        haplotypes = pd.DataFrame(data)
        expected_theta = 1.5  # Example expected value (manually calculated)
        theta_result = calc_watterson(haplotypes)
        self.assertAlmostEqual(theta_result, expected_theta, places=3)
    
    def test_calc_watterson_empty(self):
        """
        Test calc_watterson with an empty DataFrame (edge case).
        """
        empty_haplotypes = pd.DataFrame()
        self.assertEqual(calc_watterson(empty_haplotypes), 0)
    
    # You can keep adding more tests for other functions...

if __name__ == '__main__':
    unittest.main()
