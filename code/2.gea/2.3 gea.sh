############################   lfmm #########################################
library(gradientForest)
library(data.table)
gfdata <- read.csv(file = "env.csv", header = TRUE, sep=",", row.names=1)
envGF <- gfdata[,c(2:17)]
all_snp <- fread("core.lfmm",header=F)
all_snp=t(all_snp)
colnames(all_snp) <- as.character(unlist(all_snp[1, ]))
all_snp <- all_snp[-1, ]
all_snp <- all_snp[,-1]
all_snp <- all_snp[, colSums(is.na(all_snp)) == 0]
all_snp_df <- as.data.frame(all_snp)
maxLevel <- log2(0.368*nrow(all_snp_df)/2)
a <- cbind(envGF, all_snp_df)
all_gfmod <- gradientForest(a, predictor.vars = colnames(envGF),
                            response.vars = colnames(all_snp), 
                            ntree = 500, compact = TRUE, nbin = 1001,
                            maxLevel = maxLevel, trace = TRUE, 
                            corr.threshold = 0.5)
pdf(file = "sort2.pdf", width = 10, height = 15)
plot(all_gfmod, plot.type = "Overall.Importance", col = c(rep("grey", 15),
cor),las = 2, cex.axis = 0.8)
dev.off()

############################   lfmm #########################################
library(lfmm)
library(data.table)
Ya <- fread("core.lfmm", header = FALSE)
Xa <- fread("core.env", header = F, sep = ",")
enva <- Xa[, 3:12]
Ya_matrix <- as.matrix(Ya)
Xa_matrix <- as.matrix(enva)
mod.lfmma <- lfmm_ridge(Y = Ya_matrix, X = Xa_matrix, K = 1)
pv <- lfmm_test(Y = Ya_matrix, X = Xa_matrix, lfmm = mod.lfmma, calibrate = "gif")
pvalues <- pv$calibrated.pvalue
fwrite(data.table(pvalues), file ="core.pvalue")
library(qvalue)
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x,fdr.level = 0.05)
  q <- q$qvalues
  return(q)
}
p <- fread("core.pvalue")
qvalues <- apply(p, 2, my_qvalue)
fwrite(qvalues,"core.q")
############## This code is adapted from RDA analysis in *Populus koreana*#################
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
gen <- fread("core.lfmm", header = F)
LOCI <- fread("snp.id", header = F)
env <- read.table(paste0(i, ".env"), head = F, sep = ",")
env$V2 <- as.character(env$V2)
rownames(gen) <- as.character(env$V1)
colnames(gen) <- as.character(LOCI$V1)
pred <- env[, c(3:12)]
ade.rda <- rda(gen ~ ., data = pred, scale = T)
load.rda <- scores(ade.rda, choices = c(1:6), display = "species")
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}
cand1 <- outliers(load.rda[, 1], 3)
cand2 <- outliers(load.rda[, 2], 3)
cand3 <- outliers(load.rda[, 3], 3)
cand4 <- outliers(load.rda[, 4], 3)
cand5 <- outliers(load.rda[, 5], 3)
cand6 <- outliers(load.rda[, 6], 3)
cand1 <- cbind.data.frame(rep(1, times = length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, times = length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, times = length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4, times = length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5, times = length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6, times = length(cand6)), names(cand6), unname(cand6))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4) <- colnames(cand5) <- colnames(cand6) <- c("axis", "snp", "loading")
cand <- rbind(cand1, cand2, cand3, cand4, cand5, cand6)
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow = nrow(cand), ncol = 10)
gen <- as.data.frame(gen)
for (i in 1:length(cand$snp)) {
 if (i %% 10000 == 0) {
    print(i)
  }
  nam <- cand[i, 2]
  snp.gen <- gen[, nam]
  foo[i, ] <- apply(pred, 2, function(x) cor(x, snp.gen))
}
cand <- cbind.data.frame(cand, foo)
cand <- cand[!duplicated(cand$snp), ]
for (i in 1:length(cand$snp)) {
  if (i %% 10000 == 0) {
    print(i)
  }
  bar <- cand[i, ]
  cand[i, 14] <- names(which.max(abs(bar[4:13])))  
  cand[i, 15] <- max(abs(bar[4:13]))               
}
colnames(cand)[14] <- "predictor"
colnames(cand)[15] <- "correlation"
write.table(cand, file = "core_cand.txt", quote = F)
}
