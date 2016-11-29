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
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
      
class prediction:

      def RF(X, Y, cv):
            """ Predict trait using Random Forest """
            
            numberTrees = 100
            mse = []
            Y_TRUE = []
            Y_PRED = []

            for train_index, test_index in kf.split(X):
              #print("TRAIN:", train_index, "TEST:", test_index)
              X_train, X_test = X[train_index], X[test_index]
              y_train, y_test = Y[train_index], Y[test_index]
              
              model = RandomForestRegressor(n_estimators=numberTrees, max_features="sqrt")
              model.fit(X_train, y_train)

              y_pred = model.predict(X_test)

              Y_TRUE = np.append(Y_TRUE, y_test)
              Y_PRED = np.append(Y_PRED, y_pred)
              mse = np.append(mse, (mean_squared_error(y_test, y_pred)))

            Y_TRUE = np.array(Y_TRUE)
            Y_PRED = np.array(Y_PRED)
            mse = np.array(mse)
            print("KF MSE: %0.2f (+/- %0.2f)" % (mse.mean(), mse.std()*2))

            plt.scatter(Y_TRUE, Y_PRED)
            plt.xlabel('True Trait Value')
            plt.ylabel('Predicted Trait Value')
            plt.show()




            
      def NN(X, Y, cv):
            """ Predict trait using Neural Networks """
            Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=0.9, random_state=0)
            print(Xtest)   


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

      ### Data pre-processing ### 
      genotype_file = pd.read_csv(G, index_col = "Entry")
      phenotype_file = pd.read_csv(P)

      # Average pheotype values for each line
      grouped = phenotype_file.groupby(['GHID'])[T].agg([np.average]).reset_index()
      grouped = grouped.set_index("GHID")

      # Merge genotype and phenotype files by line name (GHID/Entry)
      df = pd.concat([grouped, genotype_file], axis=1, join='inner')

      # Drop rows that don't have trait value
      df = df.dropna(subset=["average"], how = "any")

      # Make X & Y for machine learning
      X = df.drop('average', axis=1).values  
      Y = df.loc[:, 'average'].values

      kf = KFold(n_splits=10)   # set k-fold number

      
      
      ### Make predictions  ###
      
      if M == "RF" or M == "RandomForest":
            from sklearn.ensemble import RandomForestRegressor
            prediction.RF(X, Y, kf)
      
      elif M == "NN" or M == "NeuralNetworks":
            prediction.NN(X, Y, kf)
      else:
            print("Prediction method not available in this script")
          
      
       ### Score Predictions ###
          
      
