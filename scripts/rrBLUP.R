#install.packages("rrBLUP")
library(rrBLUP)
library(SuppDists)
#install.packages("BGLR")
library(BLR)



####### Default values #######
#setwd("~/Desktop/Combined_stress/ARuleMining/")
num <- 12
measure <- 'lift'
file <- "NNU_Ras_df.txt_decisiontree_50"
support <- as.numeric(0.15)
confidence <- as.numeric(0.5)
###############

args = commandArgs(TRUE)
setwd(args[1])
geno_file <- args[2]
pheno_file <- as.integer(args[3])
pheno <- args[4]
sets <- args[5]


# Set of 599 wheat lines genotyped at 1279 DArT markers from the BLR package (Pérez et al. 2010)
data(wheat)

# Convert markers from {1,0} to {-1,1} - i.e.{aa,Aa,AA} for biallelic SNPs
M <- 2*X-1

# rr-BLUP using all 1279 markers - with shrinkage = 3%
A1 <- A.mat(M,shrink=TRUE)
f<-c()
for(i in 1:10){
  test <- which(sets==i)
  yNA <- Y[,1] #grain yield in environment 1
  yNA[test] <- NA #mask yields for validation set
  data1 <- data.frame(y=yNA,gid=1:599)
  rownames(A1) <- 1:599
  ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y")
  cor <-round(cor(ans1$g[test],Y[test,1]),2)
  f <- c(f,cor)
}
f

# RKHS (or rr with a linear kernel)
D <- as.matrix(dist(M)) #Euclidean distance

RKHS <-c()
for(i in 1:3){
  test <- which(sets==i)
  yNA <- Y[,1] #grain yield in environment 1
  yNA[test] <- NA #mask yields for validation set
  data1 <- data.frame(y=yNA,gid=1:599)
  rownames(A1) <- 1:599
  ans2 <- kin.blup(data1,K=D,GAUSS=TRUE, geno="gid",pheno="y",  n.core=10)
  cor <-round(cor(ans2$g[test],Y[test,1]),2)
  RKHS <- c(RKHS,cor)
}
RKHS





############### McCouch Rice Data
rice_g <- read.csv("Rice_DivPanel/MET_crfilt_.90_outliers_removed_for_RRBlup_line_corrected.csv", header=TRUE)

rice_p09 <- read.csv("Rice_DivPanel/corrected_RYT2009DS_plotdata_by_GHID.csv")
colnames(rice_p09) <- paste("09", colnames(rice_p09), sep = "_")
rice_p10 <- read.csv("Rice_DivPanel/corrected_RYT2010DS_plotdata_by_GHID.csv")
colnames(rice_p10) <- paste("10", colnames(rice_p10), sep = "_")
rice_p11 <- read.csv("Rice_DivPanel/corrected_RYT2011DS_plotdata_by_GHID.csv")
colnames(rice_p11) <- paste("11", colnames(rice_p11), sep = "_")
rice_p12 <- read.csv("Rice_DivPanel/corrected_RYT2012DS_plotdata_by_GHID.csv")
colnames(rice_p12) <- paste("12", colnames(rice_p12), sep = "_")

rice_yld <- merge(rice_p09, rice_p10, by.x=c("09_GHID", "09_REP"), by.y=c("10_GHID", "10_REP"), all.x=FALSE, all.y=FALSE)
rice_yld <- merge(rice_yld, rice_p11, by.x=c("09_GHID", "09_REP"), by.y=c("11_GHID", "11_REP"), all.x=FALSE, all.y=FALSE)
rice_yld <- merge(rice_yld, rice_p12, by.x=c("09_GHID", "09_REP"), by.y=c("12_GHID", "12_REP"), all.x=FALSE, all.y=FALSE)



pheno <- "Plot_Yld"
setwd("~/Desktop/GS_Datasets/")
g <- read.csv("genotype_sub_1200Kb_1.csv", header = TRUE)
g <- read.csv("geno_sub_120Kb_binned_random_selection_w_maf_filter_1.csv", header = TRUE)


names(g)[1] <- "GHID"
ds <- read.csv("phenotype_RYT2012DS.csv", header = TRUE)
ws <- read.csv("phenotype_RYT2012WS.csv", header = TRUE)

ds_ag <- aggregate(.~GHID, data=ds, mean)
ws_ag <- aggregate(.~GHID, data=ws, mean)

trait_of_interest <- ds_ag[,c("GHID", pheno)]

data <- merge(trait_of_interest, g, by = "GHID")

# rr-BLUP 
A1 <- A.mat(data[3:ncol(data)],shrink=TRUE)
f<-c()
for(i in 1:10){
  test <- which(map==i)
  yNA <- data$Plot_Yld #grain yield in environment 1
  yNA[test] <- NA #mask yields for testing set
  data1 <- data.frame(y=yNA, gid = 1:(nrow(data))) 
  rownames(A1) <- 1:nrow(data)
  ans1 <- kin.blup(data1,K=A1,geno="gid",pheno="y")
  cor <-round(cor(ans1$pred[test],data$Plot_Yld[test]),2)
  training_cor <-round(cor(ans1$pred[-test],data$Plot_Yld[-test]),2)
  f <- c(f,cor)
}
mean(f)
mean(training_cor)


# RKHS

D <- as.matrix(dist(data[3:ncol(data)])) #Euclidean distance

RKHS <-c()
for(i in 1:10){
  test <- which(map==i)
  yNA <- data$Plot_Yld #grain yield in environment 1
  yNA[test] <- NA #mask yields for validation set
  data1 <- data.frame(y=yNA,gid=1:(nrow(data)))
  rownames(A1) <- 1:nrow(data)
  ans2 <- kin.blup(data1,K=D,GAUSS=TRUE, geno="gid",pheno="y",  n.core=10)
  cor <-round(cor(ans1$pred[test],data$Plot_Yld[test]),2)
  RKHS <- c(RKHS,cor)
}
mean(RKHS)

