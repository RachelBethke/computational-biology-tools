#imports
import numpy as np
import pandas as pd
from itertools import combinations
import random
#from scipy.stats import tajima_d


# In[28]:
# Input matrices 
def load_haplotype_data(file_path):
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

# Initialize a list to store each haplotype matrix
haplotype_blocks = []

# Need to do the following for all 4 blocks selected from the 4 files.
#hapmx = pd.read_csv('ConGensData1block1.txt', sep = ',', header=None)
hapmx1 = load_haplotype_data('Data/ConsGenData1.txt')
hapmx2 = load_haplotype_data('Data/ConsGenData2.txt')
hapmx3 = load_haplotype_data('Data/ConsGenData3.txt')
hapmx4 = load_haplotype_data('Data/ConsGenData4.txt')

haplotype_blocks.extend([hapmx1, hapmx2, hapmx3, hapmx4])

print(haplotype_blocks[0]) #display first matrix for reference


# In[29]:
# Calculate nucelotide diversity (pi)
def calc_pi(haplotypes):
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
    
    #formula
    pi = pi_total / L  # avg nucleotide diversity
    return pi

for i, haplotype_block in enumerate(haplotype_blocks):
    result_pi = calc_pi(haplotype_block)
    #edited to look at each block
    print(f"Nucleotide diversity (pi) for block {i+1}: {result_pi}")


# In[31]:
# Calculate the number of segregating sites (S)
def count_segregating_sites(haplotypes):
    n = len(haplotypes.axes[1])
    S = 0
    for i in range(0,n):
        if len(haplotypes[i].unique()) == 1: 
            S = S
        else:
            S = S + 1
    return S

for i, hapmx in enumerate(haplotype_blocks):
    result_S = count_segregating_sites(hapmx)
    print(f"Number of segregating sites (S) for block {i+1}: {result_S} out of {haplotype_block.shape[1]} sites.")


# In[32]:
# Calculate Wattersons theta
def calc_watterson(haplotypes):
    n = haplotypes.shape[0]  # number of haplotypes (rows)
    L = haplotypes.shape[1]  # number of loci (columns)
    # S = segregating sites
    S = count_segregating_sites(haplotypes)
    # a = sum(1/i)
    a = 0
    for i in range(1,n):
        a += (1/i)
    # pi = S / a ...
    if a != 0:
        theta = S / a
    else:
        theta = 0
    return theta

for i, hapmx in enumerate(haplotype_blocks):
    result_watt = calc_watterson(hapmx)
    print(f"Wattersons theta for block {i+1}: {result_watt}")


# In[36]:
# Calculate Tajima's D
def tajimas_d(haplotypes):
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

# made it a for loop to show all the blocks
for i, hapmx in enumerate(haplotype_blocks):
    result_tajimas_d = tajimas_d(hapmx)
    print(f"Tajima's D for block {i+1}: {result_tajimas_d}")