#########################################
library(vegan)
library(geosphere)
env=read.csv("env.csv",head=T)
fst=read.csv("fst.csv",head=T)
geo=read.csv("geo.csv",head=T)
fst_matrix <- as.dist(fst)  
geo_dist <- dist(geo, method="euclidean") 
env_dist <- dist(env, method="euclidean")  # 计算欧式距离
##### IBD
mantel(fst_matrix,geo_dist,method="pearson",permutations=999)
##### IBE
mantel(fst,env_dist,method="pearson",permutations=999)
