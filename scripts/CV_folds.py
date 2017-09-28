""" 

Calls R script to do clustering by population structure. Then creates cv folds.

export PATH=/Users/cbazodi/anaconda3/bin:$PATH

"""
import sys, os
import random
import pandas as pd
import numpy as np
from collections import defaultdict


def cluster(geno, wd):
  """ Run fpc.R script to determine ideal k (pamk) and then cluster (pam) """

  R=('R --vanilla --slave --args '+wd+' '+geno+' < /Users/cbazodi/Desktop/Genomic_Selection/scripts/fpc.R')
  os.system(R)
  print("Clustering using fpc.R is complete")
  return wd+geno+"_pams.csv"



def ps_cv_folds(df, pams, cv_folds, i):
  """ Reads in the clustering .csv and builds cv folds with close to even numbers in each """
  
  # Read in clustering csv and make dict with PAM: [Entry1, Entry2, etc.]
  cluster_dict = defaultdict(list)
  num_lines = 0
  with open(pams, 'r') as pam_file:
    for l in pam_file:
      num_lines += 1
      line, cluster = l.strip().split('\t')
      cluster_dict[int(cluster)].append(line)

  ave_size = num_lines/cv_folds

  # Make list of tuples [(PAM#, size), (PAM#, size), etc.] and sort by size
  PAM_sizes = []
  for key in cluster_dict:
    PAM_sizes.append([key,len(cluster_dict[key])])
  sort_by_size = sorted(PAM_sizes, key=lambda tup: -tup[1])
  
  cv_assignment = defaultdict(list)
  for first in range(0, cv_folds):
    cv_assignment[first+1].append(sort_by_size[0])
    del sort_by_size[0]


  # Make a dict of tuples {1: [(PAM#, size), (PAM#, size), etc.]} so that total size ~ave_size
  for fold in range(1, cv_folds+1):
    size = cv_assignment[fold][0][1]
    try_these = sort_by_size[:]
    #try:
    
    while size < ave_size and len(try_these) > 0:
      test = random.choice(try_these)
      try_these.remove(test)
      test_size = int(test[1])
      if size + test_size <= ave_size*1.02:
        cv_assignment[fold].append(test)
        size = size + test_size
        sort_by_size.remove(test)
      
  
  # If a pam cluster doesn't fit into a fold, add it to the last fold (that one should be the smallest!)
  if len(sort_by_size) > 0:
    for x in sort_by_size:
      cv_assignment[10].append(x)
  
  # Make the final key with {1: [PAM, PAM, PAM], 2: [PAM, PAM], etc}
  key = {}
  for f in cv_assignment:
    f_all = cv_assignment[f]
    for x in f_all:
      key[x[0]] = f
  
  col_name = 'cv_' + str(i)
  df[col_name] = df['pam_cluster']
  df[col_name] = df[col_name].replace(to_replace = key)

  save_name = pams[:-4]+"_CVs.csv"
  df.to_csv(save_name, sep="\t", float_format='i')

  

######## Start Script ########## 

if len(sys.argv) <= 1:
  print(__doc__)
  exit()

num = 10
cv_folds = 10

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-geno":            # data frame file (Col 1 = ID, Col2 = Class, Col3-... = Features)
    geno = sys.argv[i+1]
  if sys.argv[i] == "-wd":            # working directory
    wd = sys.argv[i+1]
  if sys.argv[i] == '-ps':            # T/F. Consider population structure in making cv folds
    ps = sys.argv[i+1]
  if sys.argv[i] == "-n":             # number of cv fold sets to make (default = 1)
    num = sys.argv[i+1]



if ps == 'T' or ps == 'True':
  # Run R script to cluster lines using pam (k calucated using pamk)
  pams = cluster(geno, wd)

  # Read in the pam output file as a pd dataframe - will use this to build the cv fold lists
  df = pd.read_csv(pams, sep = '\t', header=None, index_col = 0, names = ['Line', 'pam_cluster'])

  # Create num of cv folds (usually 10)
  for i in range(1,int(num)+1):
    ps_cv_folds(df, pams, cv_folds, i)



elif ps == 'F' or ps == 'False':
  from sklearn.model_selection import ShuffleSplit
  
  # Read in the genotype file to use the entries as the skeleton for the cv fold lists
  df = pd.read_csv(geno, header=0, usecols = ['Line'])

  number = len(df['Line'].values)
  folds_list = [] 
  while len(folds_list) < number:
    folds_list.extend(['10','9','8','7','6','5','4','3','2','1'])
  while len(folds_list) > number:
    del folds_list[0]

  for i in range(1,int(num)+1):
    random.shuffle(folds_list)
    col_name = 'cv_' + str(i)
    df[col_name] = folds_list

  save_name = geno[:-4]+"_CVs.csv"
  df.to_csv(save_name, sep="\t", float_format='i', index=False)

else: 
  print("Specifiy if cv folds should be based on population structure (-ps T) or not (-ps F)")

print("done :)")

