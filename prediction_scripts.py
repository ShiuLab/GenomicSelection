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
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class prediction:

  def RF(X, Y, cv):
    """ Predict trait using Random Forest """
    
    numberTrees = 100
    
    Y_TRUE = []
    Y_PRED = []
    Y_TRUE_train = []
    Y_PRED_train = []

    mse = []
    r2 = []

    for train_index, test_index in kf.split(X):
      #print("TRAIN:", train_index, "TEST:", test_index)
      X_train, X_test = X[train_index], X[test_index]
      y_train, y_test = Y[train_index], Y[test_index]
      
      # Build and fit the RF Regressor 
      model = RandomForestRegressor(n_estimators=numberTrees, min_samples_leaf = 5, min_samples_split=20, max_depth = 10, max_features="sqrt", random_state=42)  # , 

      model.fit(X_train, y_train)
      
      # Make predictions (test and training sets)
      y_pred = model.predict(X_test)
      y_pred_train = model.predict(X_train)
      

      Y_TRUE = np.append(Y_TRUE, y_test)
      Y_PRED = np.append(Y_PRED, y_pred)
      Y_TRUE_train = np.append(Y_TRUE_train, y_train)
      Y_PRED_train = np.append(Y_PRED_train, y_pred_train)
      
      mse = np.append(mse, (mean_squared_error(y_test, y_pred)))
      r2 = np.append(r2, (r2_score(y_test, y_pred)))

    return Y_TRUE, Y_PRED, Y_TRUE_train, Y_PRED_train, mse, r2


        
  def NN(X, Y, cv):
    """ Predict trait using Neural Networks """
    print(cv)

  def FeatSel(X, Y):
    """Feature selection using DecisionTree on the whole dataframe
    Feature importance from the Random Forest Classifier is the Gini importance
    (i.e. the normalized total reduction of the criterion for the decendent nodes
      compared to the parent node brought by that feature across all trees.)
    """
    from math import sqrt
    from sklearn.ensemble import RandomForestRegressor

    n = 100   # Number of features to keep

    FS_forest = RandomForestRegressor(criterion='mse', max_features= "sqrt", n_estimators=500, n_jobs=8)
    
    #Train the model & derive importance scores
    FS_forest = FS_forest.fit(X, Y)
    importances = FS_forest.feature_importances_

    # Sort importance scores and keep top n
    print(X)
    print(X.columns)
    print(X.values)
    feat_names = list(X.columns.values)
    temp_imp = pd.DataFrame(importances, columns = ["imp"], index=feat_names) 
    indices = np.argsort(importances)[::-1]
    indices_keep = indices[0:n]
    fixed_index = []

    good = [X.columns[i] for i in indices_keep]

    X = X.loc[:,good]
    print("Features selected using DecisionTree feature selection: %s" % str(good))
    return(X)


if __name__ == "__main__":  

  FS = "False"

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
    if sys.argv[i] == "-FS":            # Set to "True" if you want to perform feature selection using decision tres
      FS = sys.argv[i+1]    

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

  kf = KFold(n_splits=10, random_state = 42)   # set k-fold number

  
  if FS == "True" or FS == "T" or FS == "true": 
    X = prediction.FeatSel(X,Y)
  
  ### Make predictions  ###
  
  if M == "RF" or M == "RandomForest":
    from sklearn.ensemble import RandomForestRegressor
    Y_TRUE, Y_PRED, Y_TRUE_train, Y_PRED_train, mse, r2 = prediction.RF(X, Y, kf)

  elif M == "NN" or M == "NeuralNetworks":
    prediction.NN(X, Y, kf)
  else:
    print("Prediction method not available in this script")
      
  
  ### Score Predictions ###
      
  #for i in range(1, len(Y_TRUE)):
  #  print(Y_TRUE[i], Y_PRED[i])

  print("MSE: %0.2f (+/- %0.2f)" % (mse.mean(), mse.std()*2))
  print("R^2: %0.2f (+/- %0.2f)" % (r2.mean(), r2.std()*2))
  print("Accuracy testing (PCC): %0.2f" % (np.corrcoef(Y_TRUE, Y_PRED)[0,1]))
  print("Accuracy training (PCC): %0.2f" % (np.corrcoef(Y_TRUE_train, Y_PRED_train)[0,1]))
  
  #fig, ax = plt.subplots()

  #ax.scatter(Y_TRUE, Y_PRED)
  #ax.xlabel('True Trait Value')
  #ax.ylabel('Predicted Trait Value')

  #plt.show()

