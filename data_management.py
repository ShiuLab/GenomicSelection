"""
PURPOSE:
Run prediction models on GS data 

If running in HPC you must set path to Miniconda in HPC:  export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH


INPUT:
  -geno       Genotype file
  -pheno      Phenotype file
  -trait      Name of column in phenotype file with trait to predict 
  -m          Prediction model to run (RF, NN, etc.)

OUTPUT:
"""

import sys
from sklearn.cross_validation import train_test_split
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
      
class prediction:
      
      def RF(self, X, Y):
            Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=10000, random_state=0)
            print(Xtrain)
            
      
if __name__ == "__main__":      
      if len(sys.argv) <= 1:
            print(__doc__)
            exit()
      for i in range (1,len(sys.argv),2):
            if sys.argv[i] == "-geno":          # Path to Genotype File
              G = sys.argv[i+1]
            if sys.argv[i] == "-pheno":         # Path to Phenotype File
              P = sys.argv[i+1] 
            if sys.argv[i] == "-trait":         # Column name for trait to use in Phenotype File
              T = sys.argv[i+1]
            if sys.argv[i] == "-m":             # Prediction method to use
              M = sys.argv[i+1]      

      ### Data pre-processing
      
      genotype_file = pd.read_csv(G, index_col = "Entry")
      phenotype_file = pd.read_csv(P)

      # Average pheotype values for each line
      grouped = phenotype_file.groupby(['GHID'])[T].agg([np.average]).reset_index()
      grouped = grouped.set_index("GHID")

      # Merge genotype and phenotype files by line name (GHID/Entry)
      df = pd.concat([grouped, genotype_file], axis=1, join='inner')

      #print(df.head(3))     

      # Make X & Y for machine learning
      X = df.drop('average', axis=1).values  
      Y = df.loc[:, 'average'].values

      if m == "RF" or m == "RandomForest":
            prediction.RF(X,Y)
             
      else:
            print("Prediction method not available in this script")
      
