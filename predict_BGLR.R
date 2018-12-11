#####################################
# Make phenotype predictions using BGLR
#
# Arguments: [1]   X_file
#            [2]   Y_file
#            [3]   Features to keep (or "all")
#            [4]   trait (col_name or "all")
#            [5]   BGLR_model (BL, BRR, BayesA, BayesB)
#            [6]   hold out file
#            [7]   save name (also used as identifier in BGLR_RESULTS.csv file)
#            [8]   save directory
#
# Example run: 
# 
#####################################
library(BGLR)
library(data.table)


# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args) < 8) {
  stop("Need 8 arguments: X_file Y_file Feats trait BGLR_model holdout_file save_name save directory]", call.=FALSE)
} 


X_file <- args[1]
Y_file <- args[2]
feat_file <- args[3]
trait <- args[4]
BGLR_model <- args[5]
ho_file <- args[6]
save_name <- args[7]
save_dir <- args[8]
df0 <- 5


## load the phenotypes and PCs
print('Loading data...')
X <- fread(X_file)
X <- as.data.frame(X)
row.names(X) <- X$ID
X$ID <- NULL

Y <- fread(Y_file)
Y <- as.data.frame(Y)
row.names(Y) <- Y$ID
Y$ID <- NULL

ho <- scan(ho_file, what='character')

# Make sure Y is in the same order as X:
Y <- Y[match(rownames(X), rownames(Y)),]
feat_method <- 'none'

# Subset X if feat_file is not all
if (feat_file != 'all'){
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character')
  X <- X[FEAT]
  feat_method <- tail(unlist(strsplit(feat_file, '/')), n=1)
  feat_method <- unlist(strsplit(feat_method, '_'))[3]
}


feat_num <- dim(X)[2]
ho_name <- tail(unlist(strsplit(ho_file, '/')), n=1)


if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}


for(i in 1:length(Y)){
  print('Building Model...')
  y=data.frame(Y[, names(Y)[i]])
  row.names(y) <- row.names(Y)
  names(y) <- c('y')
  
  # Mask yields for holdout set
  yNA <- y
  yNA[ho,] <- NA 
  yNA <- unlist(yNA)

  ETA=list(list(X=X,model=BGLR_model)) 
  fm=BGLR(y=yNA,ETA=ETA,verbose=FALSE,nIter=12000,burnIn=2000)
  y$yhat <- fm$yHat
  test_res <- y[ho,]
  coef <- fm$ETA[[1]]$b

  # Calculate scoring metrics
  pcc <- cor(test_res$yhat, test_res$y)
  mse <- mean((test_res$y - test_res$yhat)^2)
  
  # Save predicted Y
  yhat_out <- paste(BGLR_model, save_name, names(Y)[i], 'yhat.csv', sep='_')
  job_ID <- paste(BGLR_model, save_name, names(Y)[i], ho_name, sep='_')
  colnames(y) <- c('y', job_ID)
  yhat_t <- t(y)
  yhat_t <- yhat_t[-c(1), ,drop=FALSE]
  write.table(yhat_t, yhat_out, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(yhat_out))
  
  # Save coefficients  results files
  coef_out <- paste(BGLR_model, save_name, names(Y)[i], 'coef.csv', sep='_')
  coef_df <- as.data.frame(coef)
  names(coef_df) <- job_ID
  coef_df <- t(coef_df)
  write.table(coef_df, coef_out, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(coef_out))
  
  # Add to BGLR_RESULTS.csv file
  time.taken <- difftime(Sys.time(), start.time, units='sec')
  res <- data.frame(BGLR_model, X_file, save_name, names(Y)[i], ho_name, feat_method, feat_num, pcc, pcc^2, mse, time.taken)
  colnames(res) <- c('model', 'x_file','tag', 'y', 'holdout_set','feat_method','feat_num', 'PCC', 'r2','mse', 'run_time')
  write.table(res, 'BGLR_RESULTS.csv', sep=',', append=T, row.names=F, quote=F, col.names=!file.exists('BGLR_RESULTS.csv'))  
  

}

unlink('*.dat')

print('Done!')
