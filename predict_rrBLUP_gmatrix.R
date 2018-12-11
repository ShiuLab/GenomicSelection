#####################################
# Make phenotype predictions using rrBLUP
#
# Arguments: [1]   X_file
#            [2]   Y_file
#            [3]   Features to keep
#            [4]   trait (col_name or all)
#            [5]   Hold out set
#            [6]   Save name
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
if (length(args) < 6) {
  stop("Need 6 arguments: X_file Y_file Feats trait CVFs_file cvJobNum save_name [optional: save directory]", call.=FALSE)
} else if (length(args) < 7) {
  # default output file
  args[7] <- "/mnt/home/azodichr/03_GenomicSelection/"
}
#setwd('/Volumes/azodichr/03_GenomicSelection/spruce_Dial_beaulieu/')
#X_file <- '01_Data/geno.csv'
#Y_file <- '01_Data/pheno.csv'
#ho_file <- '11_FeatureSel/holdout_set/holdout2.txt'
#feat_file <- '11_FeatureSel/01_FS_bayesA/spruce_HT_bayesa_2_10'
#trait <- 'HT'
#jobNum <- 2
#save_dir <- setwd('/Volumes/azodichr/03_GenomicSelection/sorgh_DP_Fernan/11_FeatureSel/03_rrBLUP/')

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
ho_num <- tail(unlist(strsplit(ho_file, '/')), n=1)
ho_num <- gsub('holdout', '', ho_num)
ho_num <- gsub('.txt', '', ho_num)
ho_num <- gsub('.csv', '', ho_num)
ho_num <- gsub('_', '', ho_num)

# Make the relationship matrix from the markers
M=tcrossprod(scale(X))  # centered and scaled XX'
M=M/mean(diag(M))
rownames(M) <- 1:nrow(X)

if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}


# Make output directory
setwd(save_dir)


for(i in 1:length(Y)){
  print('Building Model...')
  y=data.frame(Y[, names(Y)[i]])
  row.names(y) <- row.names(Y)
  names(y) <- c('y')
  
  # Mask yields for holdout set
  yNA <- y
  yNA[ho,] <- NA 
  
  # Set up dataframe with traits and genotype labels (same order as in M) 
  df <- data.frame(y=yNA,gid=1:nrow(X)) 
  names(df) <- c('y', 'gid')
  
  # Build rrBLUP model and save yhat for the masked values
  rrblup <- kin.blup(df,K=M,geno="gid",pheno='y') 
  y$yhat<- rrblup$g
  
  holdout_pred <- y[ho,]
  pcc <- cor(holdout_pred$yhat, holdout_pred$y)
  mse <- mean((holdout_pred$y - holdout_pred$yhat)^2)
  

  # Predicted Y
  yhat_out <- paste('rrBLUP', save_name, names(Y)[i], 'yhat.csv', sep='_')
  job_ID <- paste('rrBLUP', save_name, names(Y)[i], ho_num, sep='_')
  colnames(y) <- c('y', job_ID)
  yhat_t <- t(y)
  yhat_t <- yhat_t[-c(1), ,drop=FALSE]
  write.table(yhat_t, yhat_out, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(save_name))
  
  # save  results files
  time.taken <- difftime(Sys.time(), start.time, units='sec')
  
  # Composite accuracy.txt output
  res <- data.frame('rrBLUPgmatrix', X_file, save_name, names(Y)[i], ho_name, feat_method, feat_num, pcc, pcc^2, mse, time.taken)
  colnames(res) <- c('model', 'x_file','tag', 'y', 'holdout_set','feat_method','feat_num', 'PCC', 'r2','mse', 'run_time')
  write.table(res, 'rrBLUPgmatrix_RESULTS.csv', sep=',', append=T, row.names=F, quote=F, col.names=!file.exists('rrBLUP_RESULTS.csv'))  
}

unlink('*.dat')
  


print('Complete')
