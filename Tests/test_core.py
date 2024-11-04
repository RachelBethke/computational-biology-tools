import unittest
import pandas as pd
import numpy as np
from src.genetics.core import calc_pi, calc_watterson, tajimas_d 

class TestGenomicsCalculations(unittest.TestCase):
    
    def test_calc_pi_basic(self):
        """
        Test the calc_pi function with a basic haplotype block.
        """
        # example haplotype block (4 individuals, 3 loci)
        data = {
            0: [0, 1, 0, 1],
            1: [1, 0, 1, 0],
            2: [0, 1, 1, 0]
        }
        haplotypes = pd.DataFrame(data)
        expected_pi = (0.6667 + 0.6667 + 0.6667) / 3 # Corrected from previously incorrect calculations 
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
        # p = 2, h_k = 2*2*(4-2)/(4*3) = 0.6667 (corrected calculation)
        expected_pi = 0.6667
        pi_result = calc_pi(haplotypes)
        self.assertAlmostEqual(pi_result, expected_pi, places=3)
    
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
        expected_theta = 1.636  # Adjusted for precision
        theta_result = calc_watterson(haplotypes)
        self.assertAlmostEqual(theta_result, expected_theta, places=3)
    
    def test_calc_watterson_empty(self):
        """
        Test calc_watterson with an empty DataFrame (edge case).
        """
        empty_haplotypes = pd.DataFrame()
        self.assertEqual(calc_watterson(empty_haplotypes), 0)
    
    def test_calc_watterson_one_individual(self):
        """
        Test calc_watterson with only one individual (n=1).
        """
        data = {
            0: [0],
            1: [1],
            2: [1]
        }
        haplotypes = pd.DataFrame(data)
        expected_theta = 0  # 0 bc there is only one individual
        theta_result = calc_watterson(haplotypes)
        self.assertEqual(theta_result, expected_theta)

    def test_calc_watterson_identical_haplotypes(self):
        """
        Test calc_watterson with identical haplotypes (no segregating sites).
        """
        data = {
            0: [1, 1, 1, 1],
            1: [1, 1, 1, 1],
            2: [1, 1, 1, 1]
        }
        haplotypes = pd.DataFrame(data)
        expected_theta = 0 #no segregating sites makes theta 0
        theta_result = calc_watterson(haplotypes)
        self.assertEqual(theta_result, expected_theta)


    def test_calc_pi_varied_allele_frequencies(self):
        """
        Test calc_pi with a more complex haplotype block where the allele frequencies vary.
        """
        data = {
            0: [0, 1, 0, 1],
            1: [1, 1, 0, 0],
            2: [1, 1, 1, 0],
            3: [0, 0, 1, 0]
        }
        haplotypes = pd.DataFrame(data)
        # corrected expected value based on the actual data
        expected_pi = (0.6667 + 0.6667 + 0.5 + 0.5) / 4
        pi_result = calc_pi(haplotypes)
        self.assertAlmostEqual(pi_result, expected_pi, places=3)

    def test_tajimas_d_no_segregating_sites(self):
        """
        Test Tajima's D when there are no segregating sites (S = 0).
        """
        data = {
            0: [0, 0, 0],
            1: [0, 0, 0],
            2: [0, 0, 0]
        }
        haplotypes = pd.DataFrame(data)
        expected_tajimas_d = 0 #no segregating sites makes it 0
    
        tajimas_d_result = tajimas_d(haplotypes)
        self.assertEqual(tajimas_d_result, expected_tajimas_d)

    def test_tajimas_d_small_sample(self):
        """
        Test Tajima's D with a small sample size (n=2).
        """
        data = {
            0: [0, 1],
            1: [1, 0]
        }
        haplotypes = pd.DataFrame(data)
        expected_tajimas_d = 0 # b/c of small sample size it rounds down
        
        tajimas_d_result = tajimas_d(haplotypes)
        self.assertAlmostEqual(tajimas_d_result, expected_tajimas_d, places=3)

    def test_large_haplotype_data(self):
        """
        Test the calculations with a large dataset (e.g., 1000 individuals, 100 loci).
        """
        np.random.seed(0)  # for consistancy
        large_haplotypes = pd.DataFrame(np.random.randint(0, 2, size=(1000, 100)))
    
        # just in case
        pi_result = calc_pi(large_haplotypes)
        theta_result = calc_watterson(large_haplotypes)
        tajimas_d_result = tajimas_d(large_haplotypes)
        
        self.assertIsInstance(pi_result, float)
        self.assertIsInstance(theta_result, float)
        self.assertIsInstance(tajimas_d_result, float)

if __name__ == '__main__':
    unittest.main()
