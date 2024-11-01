#imports
import numpy as np
import pandas as pd
from itertools import combinations
import random
#from scipy.stats import tajima_d

# In[28]:
def load_haplotype_data(file_path):
    """
    Can be used to load haplotype files and parse haplo data to be analyzed.
    
    Args:
        file_path: A file containing haplotype data in whitespace between blocks and commas seperating numerical values.
    
    Returns: pd.DataFrame: Parsed haplotype data
    """
    with open(file_path, 'r') as file:    # 'r' = read mode
        blocks = file.read().split('\n\n')   #split file in to blocks from newlines
        block = random.choice(blocks)
        processed_lines = []
        # from each file, select one block at random,
        for line in block.strip().split('\n'):
            stripped_line = line.strip()
            if stripped_line.replace(",", "").isdigit():
                processed_lines.append(stripped_line.split(','))
        return pd.DataFrame(processed_lines).astype(int) 
    # each block consists of haplotypes from 50 individuals (rows)

# In[29]:
def calc_pi(haplotypes):
    """
    Calculates nucleotide diversity (pi) for the given haplotype data.
    
    Args:
        haplotypes (pd.DataFrame): A DataFrame containing haplotype data (rows represent individuals, columns represent loci).
    
    Returns: float: The calculated nucleotide diversity (pi) value.
    """
    # variables
    n = haplotypes.shape[0]  # number of haplotypes (rows)
    L = haplotypes.shape[1]  # number of loci (columns)
    # calculations
    if L == 0:
        return 0
    pi_total = 0
    for k in range(L):
        p = np.sum(haplotypes.iloc[:, k] == 1)  # number of derived alleles (1's) at loci k
        h_k = 2 * p * (n - p) / (n * (n - 1))
        pi_total += h_k
        #print(f"Locus {k}: p = {p}, h_k = {h_k}")
    
    #formula
    pi = pi_total / L  # avg nucleotide diversity
    return pi


# In[31]:
def count_segregating_sites(haplotypes):
    """
    Can be used to calculate the number of segregating sites.
    
    Args:
        haplotypes (pd.DataFrame): A DataFrame containing haplotype data.
        
    Returns: int: the number segregating sites (S) in the given haplotype data
    """
    n = len(haplotypes.axes[1])
    S = 0
    for i in range(0,n):
        if len(haplotypes[i].unique()) == 1: 
            S = S
        else:
            S = S + 1
    return S

# In[32]:
def calc_watterson(haplotypes):
    """
    Calculates Watterson's theta for the given haplotype data.
    
    Args:
        haplotypes (pd.DataFrame): A DataFrame containing haplotype data.
    
    Returns: float: Watterson's theta value.
    """
    n = haplotypes.shape[0]  # number of haplotypes (rows)

    if n <= 1: # edge case
        return 0 
    
    L = haplotypes.shape[1]  # number of loci (columns)
    # S = segregating sites
    S = count_segregating_sites(haplotypes)
    # a = sum(1/i)
    # a = 0
    # for i in range(1,n):
    #     a += (1/i)
    a = np.sum([1/i for i in range(1, n)])
    # pi = S / a ...
    if a != 0:
        theta = S / a
    else:
        theta = 0
    return theta



# In[36]:
def tajimas_d(haplotypes):
    """
    Calculates Tajima's D for the given haplotype data.
    
    Args:
        haplotypes (pd.DataFrame): A DataFrame containing haplotype data.
    
    Returns: float: Tajima's D value.
    """
    # variables
    n = haplotypes.shape[0]  # number of haplotypes (rows)
    pi = calc_pi(haplotypes)
    theta = calc_watterson(haplotypes)
    S = count_segregating_sites(haplotypes)

    # formulas and calculations
    a1 = np.sum([1/i for i in range(1, n)])  # Sum of 1/i from 1 to n-1
    a2 = np.sum([1/(i**2) for i in range(1, n)])  # Sum of 1/i^2 from 1 to n-1

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

        # Variance of pi - theta
    variance = (e1 * S) + (e2 * S * (S - 1))
    
    sd = np.sqrt(variance)

    # Tajima's D calculation
    if sd != 0:
        tajimas_d = (pi - theta) / sd
    else:
        tajimas_d = 0
    return tajimas_d
