# GenomicSelection





# Make CV clusters based on population structure

### 1. Make Clusters 
1. Determine the number of clusters using PAMK function in the fpc package in R (Partitioning around medoids with estimation of number of clusters). Use asw (average silhouette width - how close each point in one cluster is to points in the neighboring clusters, values range from 0-1, close to 1 is better!) for the clustering criteria. Check's 1-100 clusters.
2. Make clusters using pam (Partitioning (clustering) of the data into k clusters ``around medoids''- more robust than k-means)
3. 


