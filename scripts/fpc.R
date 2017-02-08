#install.packages("fpc")
library(fpc)
library(cluster)

args = commandArgs(TRUE)
setwd(args[1])
g <- read.csv(args[2])

#setwd("~/Desktop/Genomic_Selection/GS_Datasets/")
#g <- read.csv("geno_sub_120Kb_binned_random_selection_w_maf_filter_1.csv")
#geno_sub_120Kb_binned_random_selection_w_maf_filter_1.csv_pams.csv
# Since Spindel et al already did that for the rice dataset, just use their finding of k = 87
#k <- 87

# Remove rows with NA values
g <- na.omit(g)

# To figure out ideal number of clusters
pk<- pamk(g[,2:ncol(g)], krange=2:100, criterion="asw")
k <- pk$nc
k
pams <- pam(g, k, cluster.only=TRUE)
p <- data.frame(Entry=g$Entry,PAM=pams)
pams_sort <- sort(table(pams))

name=paste(basename(args[2]),"_pams.csv", sep="")
#name="test_pams.csv"

write.table(p, file=name, row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)


