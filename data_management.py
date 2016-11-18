import sys
from sklearn.cross_validation import train_test_split

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-geno":
        G = sys.argv[i+1]
      if sys.argv[i] == "-pheno":
        P = sys.argv[i+1] 
      if sys.argv[i] == "-trait":
        T = sys.argv[i+1]

genotype_file = pd.read_csv(G)
phenotype_file = pd.read_csv(P)

# Average pheotype values for each line
grouped = phenotype_file.groupby(['GHID'])[T].agg([np.sum, np.average]).reset_index()


X = genotype_file.values  
Y = phenotype_file


print(genotype_file.head(3))
print(phenotype_file.head(4))
print(grouped.head(3))


