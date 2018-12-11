#####################################
# Make phenotype predictions using mixed solve from rrBLUP
# Equivalent to using kin.blup on the G-matrix. See: 
# https://pdfs.semanticscholar.org/5c52/d445ea806770331b65d1ac65afaa74a619c4.pdf
#
# Arguments: [1]   X_file
#            [2]   Y_file
#            [3]   Features to keep
#            [4]   trait (col_name or all)
#            [5]   Hold out set
#            [6]   save sane (also used as identifier in rrBLUP_RESULTS.csv file)
#            [7]   optional: save directory
#
# Written by: Christina Azodi
# Original: 4.26.17
# Modified: 
#####################################
library(rrBLUP)
library(data.table)


# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args) < 7) {
  stop("Need 7 arguments: X_file Y_file Feats trait CVFs_file cvJobNum save_name save directory", call.=FALSE)
}

X_file <- args[1]
Y_file <- args[2]
feat_file <- args[3]
trait <- args[4]
ho_file <- args[5]
save_name <- args[6]
save_dir <- args[7]

## load the phenotypes and PCs
print('Loading data...')
X <- fread(X_file)
X <- as.data.frame(X)
row.names(X) <- X$ID
X$ID <- NULL

print(X[1:5,1:5])
Y <- fread(Y_file)
Y <- as.data.frame(Y)
row.names(Y) <- Y$ID
Y$ID <- NULL

ho <- scan(ho_file, what='character')
train_instances <- setdiff(row.names(X), ho) 
  
# Make sure Y is in the same order as X:
Y <- Y[match(rownames(X), rownames(Y)),]
feat_method <- 'none'

# Subset X if feat_file is not all
if (feat_file != 'all'){
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character')
  X <- X[FEAT]
  feat_method <- tail(unlist(strsplit(feat_file, '/')), n=1)
}

feat_num <- dim(X)[2]
ho_name <- tail(unlist(strsplit(ho_file, '/')), n=1)

if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

# Make output directory
setwd(save_dir)


for(i in 1:length(Y)){
  print('Building Model...')
  
  test_x <- as.matrix(X[ho,])
  test_y <- Y[names(Y)[i]][ho,]
  train_x <- as.matrix(X[train_instances,])
  train_y <- Y[names(Y)[i]][train_instances,]
  mod <- mixed.solve(train_y, Z=train_x, K=NULL, method='ML', SE=F, return.Hinv=F)
  
  # Predict test set (i.e. holdout)
  print('Predicting test set...')
  coef <- mod$u
  effect_size <- as.matrix(coef)
  yhat <- (test_x %*% effect_size)[,1]
  yhat <- yhat + mod$beta
  yhat_all <- (as.matrix(X) %*% effect_size)[,1]
  yhat_all <- yhat_all + mod$beta
  
  # Calculate scoring metrics
  print('Calculating scoring metrics...')
  pcc <- cor(yhat, test_y)
  mse <- mean((test_y - yhat)^2)

  
  # Save predicted Y (i.e. yhat)
  print('Saving output...')
  yhat_out <- paste('rrBLUP', save_name, names(Y)[i], 'yhat.csv', sep='_')
  job_ID <- paste('rrBLUP', save_name, names(Y)[i], ho_name, sep='_')
  yhat_df <- data.frame(yhat=yhat_all)
  names(yhat_df) <- job_ID
  row.names(yhat_df) <- row.names(Y)
  yhat_t <- t(yhat_df)
  write.table(yhat_t, yhat_out, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(yhat_out))
  
  # Save coefficients
  coef_out <- paste('rrBLUP', save_name, names(Y)[i], 'coef.csv', sep='_')
  coef_df <- data.frame(u=coef)
  names(coef_df) <- job_ID
  coef_df <- t(coef_df)
  write.table(coef_df, coef_out, sep=',', append=T, row.names=T, quote=F, col.names=!file.exists(coef_out))
  
  # save  results files
  time.taken <- difftime(Sys.time(), start.time, units='sec')
  
  res <- data.frame('rrBLUP', X_file, save_name, names(Y)[i], ho_name, feat_method, feat_num, pcc, pcc^2, mse, time.taken)
  colnames(res) <- c('model', 'x_file','tag', 'y', 'holdout_set','feat_method','feat_num', 'PCC', 'r2','mse', 'run_time')
  write.table(res, 'rrBLUP_RESULTS.csv', sep=',', append=T, row.names=F, quote=F, col.names=!file.exists('rrBLUP_RESULTS.csv'))  
}

unlink('*.dat')
  


print('Done!')
