############## This code is adapted from RDA analysis in *Populus koreana*#################
library(data.table)
library(gradientForest)
gfdata <- read.csv(file = "cur.csv", header = T ,row.names = 1)
envGF <- gfdata[,c(2:20)]
all_snp <- fread("../core.lfmm",header=F)
maxLevel <- log2(0.368*nrow(all_snp)/2)
a <- cbind(envGF, all_snp)
all_gfmod <- gradientForest(a, predictor.vars = colnames(envGF),
                            response.vars = colnames(all_snp), 
                            ntree = 500, compact = TRUE, nbin = 1001,
                            maxLevel = maxLevel, trace = TRUE, 
                            corr.threshold = 0.5)
pdf(file="picture/4w_predictorcumulative.pdf")
plot(all_gfmod, plot.type="C", 
imp.vars=c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"), 
     show.species=F, common.scale=T, cex.axis=0.6, cex.lab=0.7, 
     line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), 
                                  mar=c(2.5,1.0,0.1,0.5), 
                                  omi=c(0,0.3,0,0)))
dev.off()
                  
pdf(file="picture/all13_R2.pdf")
plot(all_gfmod, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)
dev.off()

pdf(file="picture/all13_splitsdensityplots.pdf")
plot(all_gfmod, plot.type="S", imp.vars=c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"), leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))
dev.off()

pdf(file="picture/all13_speciescumulativeplot.pdf")
plot(all_gfmod, plot.type="Cumulative.Importance", imp.vars=c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"), show.overall=T, legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4, cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))
dev.off()

################################################################################################################
allpos= fread("ade.txt",header=T,sep="\t")
allpos =  allpos[,c(2,3)]
greengrid  = fread("ade_19_current.csv",header=T,sep=",")
greengrid[, X := NULL]
cbinded =cbind(allpos[,c("lon","lat")],greengrid[,1:19])
all_tgrid=cbind(cbinded[,c("lon","lat")],
                predict(all_gfmod,greengrid[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]))
n<-sum(is.na(all_tgrid)) 
Trns_grid <- na.omit(all_tgrid) 
n<-sum(is.na(Trns_grid)) 
n<-sum(is.na(all_tgrid))
Trns_grid <- na.omit(all_tgrid)
n<-sum(is.na(Trns_grid))
all_PCs <- prcomp(Trns_grid[, c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")], center=TRUE, scale.=FALSE)
summary(all_PCs)
a1 <- all_PCs$x[,1]
a2 <- all_PCs$x[,2]
a3 <- all_PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255
grid <- cbinded[,c("lon","lat")]
grid$R=r
grid$G=g
grid$B=b
nvs <- dim(all_PCs$rotation)[1]
vec <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19") 
lv <- length(vec)
vind <- rownames(all_PCs$rotation) %in% vec
scal <- 60
xrng <- range(all_PCs$x[,1], all_PCs$rotation[,1]/scal)*1.1
yrng <- range(all_PCs$x[,2], all_PCs$rotation[,2]/scal)*1.1

pdf(file="picture/all13_PCplot02.pdf")
plot((all_PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
all_PCs$rotation <- all_PCs$rotation * 2
arrows(rep(0,lv), rep(0,lv), all_PCs$rotation[,1]/scal, all_PCs$rotation[,2]/scal, length = 0.1)
jit <- 0.0015
text(all_PCs$rotation[,1]/scal+jit*sign(all_PCs$rotation[,1]), all_PCs$rotation[,2]/scal+jit*sign(all_PCs$rotation[,2]), labels = vec)
dev.off()

pdf("picture/all13_Map2.pdf")
green.pred <- predict(all_gfmod, cbinded[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")  ])
plot(Trns_grid[, c("lon", "lat")], pch=15,cex=1.0,asp=1,col=rgb(r,g,b, max=255),main="SNP turnover in Q rugosa")
dev.off()

greencols=rgb(r,g,b,max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)
gradients=cbind(Trns_grid[,1:2],greencols3)
gradients$color=greencols
write.csv(gradients,file="result/all13_gradients4arcgis03.csv",row.names=F,quote=F)      
all_tgrid=cbind(cbinded[,c("lon","lat")],
                predict(all_gfmod,greengrid[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]))
                
#############################################    local    ###########################################################
ly = c("ade_19_5854160_try","ade_19_5856180_try","ade_19_1264160_try","ade_19_1266180_try")
for (i in ly) {
    print(i)
    # 读取 CSV 文件
    fut = fread(paste0(i, ".csv"), header = TRUE, sep = ",")
    dim(fut)
    # 删除第一列
    fut[, X := NULL]
    
    # 将经纬度和前6列数据合并
    fut_cbinded = cbind(allpos[, c("lon", "lat")], fut[, 1:19])
    # 基于梯度森林模型计算未来的环境梯度
    future_all = cbind(fut_cbinded[, c("lon", "lat")], 
                       predict(all_gfmod, fut_cbinded[, c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")]))
    # 计算遗传偏移
    genOffsetAll = sqrt((future_all[, 3] - all_tgrid[, 3])^2 +
                        (future_all[, 4] - all_tgrid[, 4])^2 +
                        (future_all[, 5] - all_tgrid[, 5])^2 +
                        (future_all[, 6] - all_tgrid[, 6])^2 +
                        (future_all[, 7] - all_tgrid[, 7])^2 +
                        (future_all[, 8] - all_tgrid[, 8])^2 +
                        (future_all[, 9] - all_tgrid[, 9])^2 +
                        (future_all[, 10] - all_tgrid[, 10])^2 +
                        (future_all[, 11] - all_tgrid[, 11])^2 +
                        (future_all[, 12] - all_tgrid[, 12])^2 +
                        (future_all[, 13] - all_tgrid[, 13])^2 +
                        (future_all[, 14] - all_tgrid[, 14])^2 +
                        (future_all[, 15] - all_tgrid[, 15])^2 +
                        (future_all[, 16] - all_tgrid[, 16])^2 +
                        (future_all[, 17] - all_tgrid[, 17])^2 +
                        (future_all[, 18] - all_tgrid[, 18])^2 +
                        (future_all[, 19] - all_tgrid[, 19])^2 +
                        (future_all[, 20] - all_tgrid[, 20])^2 +
                        (future_all[, 21] - all_tgrid[, 21])^2)
    # 将坐标和遗传偏移值合并
    Offset = cbind(future_all[, c("lon", "lat")], genOffsetAll)
    # 修改遗传偏移计算结果列名为“offset”
    colnames(Offset)[3] = "offset"
    # 保存遗传偏移结果到 CSV 文件
    write.csv(Offset, paste0(i, "_local.csv"), quote = FALSE, row.names = TRUE)
}

##########################################  forward #########################################
library(doParallel)
library(foreach)
library(gdm)
cl <- makeCluster(80)
registerDoParallel(cl)
predNames <- colnames(envGF)
for (file in ly) {
  print(file)
  fut = fread(paste0(file, ".csv"), header = TRUE, sep = ",")
  fut[, X := NULL]
  fut_cbinded = cbind(allpos[,c("lon","lat")], fut[,1:19])
  futClimDatGF <- data.frame(fut_cbinded[,c("lon","lat")], predict(all_gfmod, fut_cbinded[,3:21])) 
  popDatGF <- data.frame(cbinded[,c("lon","lat")], predict(all_gfmod, cbinded[,3:21])) 
  popDatGF <- split(popDatGF, seq(nrow(popDatGF)))
  predNames <- colnames(envGF)
  forwardOffsetGF <- foreach(j = 1:length(popDatGF), .packages = c("fields", "gdm", "geosphere")) %dopar% {
    onePopGF <- popDatGF[[j]]
    combinedDatGF <- futClimDatGF[,c("lon","lat")]
    combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
    coordGF <- onePopGF[,c("lon","lat")]
    minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
    minCoordsGF["dists"] <- distGeo(p1 = coordGF, p2 = minCoordsGF[,1:2])
    minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
    minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF), 1),]
    offsetGF <- combinedDatGF[which(combinedDatGF$lon == coordGF$lon & combinedDatGF$lat == coordGF$lat), "gfOffset"]
    minValGF <- minCoordsGF$gfOffset
    toGoGF <- minCoordsGF$dists
    minPtGF <- minCoordsGF[,c("lon","lat")]
    bearGF <- bearing(coordGF, minPtGF)
    outGF <- c(x1 = coordGF[[1]], y1 = coordGF[[2]], local = offsetGF, forwardOffset = minValGF, predDist = toGoGF, 
               bearing = bearGF, x2 = minPtGF[[1]], y2 = minPtGF[[2]])
    return(outGF)
  }
  forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
  write.csv(forwardOffsetGF, paste0(file, "_forwardOffsetGF.csv"), row.names = FALSE)
}
#######################################       migrant #######################################################
distances = c( 50000, 100000, 200000, 500000)
for (file in ly) {
  print(file)
  fut = fread(paste0(file, ".csv"), header = TRUE, sep = ",")
  fut[, X := NULL]
  fut_cbinded = cbind(allpos[,c("lon","lat")], fut[,1:19])
  futClimDatGF <- data.frame(fut_cbinded[,c("lon","lat")], predict(all_gfmod, fut_cbinded[,3:21])) 
  popDatGF <- data.frame(cbinded[,c("lon","lat")], predict(all_gfmod, cbinded[,3:21])) 
  popDatGF <- split(popDatGF, seq(nrow(popDatGF)))
  predNames <- colnames(envGF)
  for (dist_limit in distances) {
      print(dist_limit)
    forwardOffsetGF <- foreach(i = 1:length(popDatGF), .packages=c("fields", "gdm", "geosphere")) %dopar% {
      onePopGF <- popDatGF[[i]]
      combinedDatGF <- futClimDatGF[,c("lon","lat")]
      combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
      coordGF <- onePopGF[,c("lon","lat")]
      combinedDatGF['dists'] = distGeo(p1 = coordGF, p2 = combinedDatGF[,1:2])
      combinedDatGF <- combinedDatGF[combinedDatGF$dists < dist_limit,]
      minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
      if (nrow(minCoordsGF) > 1) {
        minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF), 1),]
      }
      offsetGF <- combinedDatGF[which(combinedDatGF$lon == coordGF$lon & combinedDatGF$lat == coordGF$lat), "gfOffset"]
      minValGF <- minCoordsGF$gfOffset
      toGoGF <- minCoordsGF$dists
      minPtGF <- minCoordsGF[,c("lon","lat")]
      bearGF <- bearing(coordGF, minPtGF)
      outGF <- c(x1 = coordGF[[1]], y1 = coordGF[[2]], local = offsetGF, forwardOffset = minValGF, predDist = toGoGF, 
                 bearing = bearGF, x2 = minPtGF[[1]], y2 = minPtGF[[2]])
      return(outGF)
    }
    forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
    write.csv(forwardOffsetGF, paste0(file, "_", dist_limit/1000, "kmforwardOffsetGF_.csv"), row.names = FALSE)
  }
}
############################################# reverse  ###################################################
for (file in ly) {
  print(file)
  fut = fread(paste0(file, ".csv"), header = TRUE, sep = ",")
  fut[, X := NULL]
  fut_cbinded = cbind(allpos[,c("lon","lat")], fut[,1:19])
  popDatGF <- data.frame(cbinded[,c("lon","lat")], predict(all_gfmod, cbinded[,3:21])) 
  futClimDatGF <- data.frame(fut_cbinded[,c("lon","lat")], predict(all_gfmod, fut_cbinded[,3:21])) 
  predNames <- colnames(popDatGF)[3:ncol(popDatGF)]
  reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","gdm","geosphere")) %dopar% {
      onePopGF <- futClimDatGF[i,]
      combinedDatGF <- popDatGF[,c("lon","lat")]
      combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], popDatGF[,predNames]))
      coordGF <- onePopGF[,c("lon","lat")]
      minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
      minCoordsGF["dists"] <- distGeo(p1 = coordGF, p2 = minCoordsGF[,1:2])
      minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
      if (nrow(minCoordsGF) > 1) {
        minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF), 1),]
      }
      offsetGF <- combinedDatGF[which(combinedDatGF$lon == coordGF$lon & combinedDatGF$lat == coordGF$lat), "gfOffset"]
      minValGF <- minCoordsGF$gfOffset
      toGoGF <- minCoordsGF$dists
      minPtGF <- minCoordsGF[,c("lon","lat")]
      bearGF <- bearing(coordGF, minPtGF)
      outGF <- c(x1 = coordGF[[1]], y1 = coordGF[[2]], local = offsetGF, reverseOffset = minValGF, predDist = toGoGF, 
                 bearing = bearGF, x2 = minPtGF[[1]], y2 = minPtGF[[2]])
      return(outGF)
  }
  reverseOffsetGF <- do.call(rbind, reverseOffsetGF)
  write.csv(reverseOffsetGF, paste0(file, "_reverse_OffsetGF.csv"), row.names = FALSE)
}
