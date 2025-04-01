library(vegan)
library(data.table)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
env <- read.csv(file = "env.csv", header = TRUE, sep=",", row.names=1)
pred=env[,c(2:17)]
gen=fread("core.lfmm",header=F)
LOCI=fread("snp.id",header=F)
rownames(gen) <- rownames(env)
colnames(gen)=as.character(LOCI$V1)
pred=env[,c(2:17)]
bio=pred[,c(1:5)]
soil=pred[,c(6:10)]
mem=pred[,c(11:13)]
str=pred[,c(14:16)]
bio2=pred[,c(6:16)]
soil2=pred[,c(1:5,11:16)]
mem2=pred[,c(1:10,14:16)]
str2=pred[,c(1:13)]
bio_soil=pred[,c(1:10)]
bio_soil2=pred[,c(11:16)]
r_bio_soil2 <- rda(gen,bio_soil,bio_soil2,scale = T)
RsquareAdj(r_bio_soil2) 
all_rda <- rda(gen ~ ., data=pred, scale=T)
RsquareAdj(all_rda) 
r_bio <- rda(gen,bio,bio2,scale = T)
RsquareAdj(r_bio) 
r_soil <- rda(gen,soil,soil2,scale = T)
RsquareAdj(r_soil) 
r_mem <- rda(gen,mem,mem2,scale = T)
RsquareAdj(r_mem) 
r_str <- rda(gen,str,str2,scale = T)
RsquareAdj(r_str) 
all_rda <- rda(gen ~ ., data=pred, scale=T)
RsquareAdj(all_rda) 
all_rda_result <- anova.cca(all_rda, permutations = 999)
all_r_bio <- anova.cca(r_bio, permutations = 999)
all_r_soil <- anova.cca(r_soil, permutations = 999)
all_r_mem <- anova.cca(r_mem, permutations = 999)
all_r_str <- anova.cca(r_str, permutations = 999)
