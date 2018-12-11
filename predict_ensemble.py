"""
Take multiple _scores.txt and _yhat.csv files and 
generate an ensemble regression accuracy and yhats

Approach:
  1. For every replicate (default 100), grab the yhats from each model 
  2. Calculate the mean yhat for that replicate
  3. Calculate the correlation(ensemble yhat ~ y) and save in accuracy.txt
  4. Save all those ensemble yhats in the SAVE file

python ensemble_pred.py SAVE TRUE_Y_FILE Y_COL_NAME DTYPE_ID FILES...

After first 4 args, can give as many _scores.txt and _yhat.csv files as desired

"""


import sys, os
from sklearn.metrics import mean_squared_error, r2_score, explained_variance_score
import pandas as pd
import numpy as np
import time

items = []
SEP = '\t'
y_name = sp = ''
reps = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

for i in range (1,len(sys.argv),2):
	if sys.argv[i].lower() == "-save":
		SAVE = sys.argv[i+1]
	if sys.argv[i].lower() == "-ho":
		ho_file_string = sys.argv[i+1]
	if sys.argv[i].lower() == "-yhat":
		yhat_file = sys.argv[i+1]
	if sys.argv[i].lower() == "-y":
		y_file = sys.argv[i+1]
	if sys.argv[i].lower() == "-sep":
		SEP = sys.argv[i+1]
	if sys.argv[i].lower() == "-y_name":
		y_name = sys.argv[i+1]
	if sys.argv[i].lower() == "-sp":
		sp = sys.argv[i+1]


# If features and class info are in separate files, merge them: 
y = pd.read_csv(y_file, sep=SEP, index_col = 0)
if y_name != '':
	y = y[[y_name]]
y_mean = y.mean(axis=0)
y_std = y.std(axis=0)
y = (y - y_mean) / y_std


# 
print('Extracting predicted values from each dataset...')
df = pd.DataFrame()
with open(yhat_file) as yhat_file_open:
	for yhat_f in yhat_file_open:
		yhat_f = yhat_f.strip()
		
		if 'bglr' in yhat_f or 'rrblup' in yhat_f:
			yhat = pd.read_csv(yhat_f, sep=SEP, index_col = 0)
			tmp_info = pd.DataFrame(yhat.index.str.split("_", expand = True))
			tmp_info = pd.DataFrame(tmp_info[0].values.tolist(), columns = ['mod', 'sp', 'y', 'hox'])
			tmp_info['ho'] = tmp_info['hox'].str.extract('(\d+)').astype(int) 
			yhat['mod'] = tmp_info['mod'].values
			yhat['ho'] = tmp_info['ho'].values.astype(str)
			df = df.append(yhat)

		elif 'ML' in yhat_f:
			yhat = pd.DataFrame()
			for r in reps:
				try:
					yhat_f_use = yhat_f.replace('XX', r)
					x = yhat_f_use.strip().split('/')[-1]
					mod = x.strip().split('_')[4]
					tmp = pd.read_csv(yhat_f_use, sep= '\t', index_col = 0)
					tmp = tmp.T
					tmp['mod'] = mod
					tmp['ho'] = r
					yhat = yhat.append(tmp.loc['Mean',:])
				except:
					print('\tError with: %s' % yhat_f_use.strip().split('/')[-1])
			df = df.append(yhat)

		elif 'mlp' in yhat_f:
			yhat = pd.read_csv(yhat_f, sep='\t', index_col = 0)
			tmp_info = pd.DataFrame(yhat.index.str.split("_", expand = True))
			tmp_info = pd.DataFrame(tmp_info[0].values.tolist(), columns = ['sp', 'y', 'hox'])
			tmp_info['ho'] = tmp_info['hox'].str.extract('(\d+)').astype(int) 
			yhat['mod'] = 'mlp'
			yhat['ho'] = tmp_info['ho'].values.astype(str)
			df = df.append(yhat)

		else:
			print('Input type not recognized: %s' % yhat_f)

			'''
			yhat = pd.read_csv(yhat_f, sep=SEP, index_col = 0)
			tmp_info = pd.DataFrame(yhat.index.str.split("_", expand = True))
			tmp_info = pd.DataFrame(tmp_info[0].values.tolist(), columns = ['mod', 'sp', 'y', 'hox'])
			tmp_info['ho'] = tmp_info['hox'].str.extract('(\d+)').astype(int) 
			yhat['mod'] = tmp_info['mod'].values
			yhat['ho'] = tmp_info['ho'].values.astype(str)
			df = df.append(yhat)
			'''
		

# Calculate performance metrics for each holdout set
print('Calculating performance metrics for each rep...')
for r in reps:
	df_r = df.loc[df['ho'] == r]
	df_r = df_r.drop_duplicates(subset='mod')
	ho_file = ho_file_string.replace('XX', r)
	with open(ho_file) as ho_open:
		ho_list = ho_open.read().splitlines()	
	
	# Pull yhat for rep
	df_r = df_r[ho_list]
	df_r_norm = df_r.T
	df_r_norm = (df_r_norm - df_r_norm.mean()) / df_r_norm.std()
	df_r_norm = df_r_norm.T
	yhat_mean = df_r_norm.mean(axis =0)

	# Pull y for rep and merge
	y_r = y.loc[ho_list,:]
	res = pd.DataFrame(yhat_mean, columns=['yhat']).join(y_r)
	res.columns = ['yhat', 'y']
	
	# Calculate performance metrics
	mse = mean_squared_error(res['y'], res['yhat'])
	evs = explained_variance_score(res['y'], res['yhat'])
	r2 = r2_score(res['y'], res['yhat'])
	cor = np.corrcoef(res['y'], res['yhat'])[0,1]
	print(mse, evs, r2, cor)

	if not os.path.isfile('RESULTS_ens.txt'):
		out2 = open('RESULTS_ens.txt', 'a')
		out2.write('ID\tSpecies\tY\tAlg\tNumEnsemble\tHO\tMSE\tEVS\tr2\tPCC\n')
		out2.close()

	out2 = open('RESULTS_ens.txt', 'a')
	out2.write('%s\t%s\t%s\tEnsemble\t%i\t%s\t%f\t%f\t%f\t%f\n' % (SAVE, sp, y_name, len(df_r.index), r, mse, evs, r2, cor))


