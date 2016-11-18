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

g = pd.read_csv(G, sep=",")
p = pd.read_csv(P)

print(g.head(3))
print(p.head(4))


