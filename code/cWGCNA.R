rm(list=ls())
library(WGCNA)
enableWGCNAThreads()

##Set working directory to code directory
load("../data_provided/Data_forConsensus_Ctx6M.rda")


## Outlier detection and removal (for each subset separately)

sdout <- 3
for (i in 1:nSets) {
  ## Remove outliers
  ##Calculate signed, weighted biweight midcorrelation
  normadj <- (0.5+0.5*bicor(t(multiExpr[[i]]$data)))^2

  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
  print(rownames(multiExpr[[i]]$data)[outliers])
  print(table(outliers))
  # multiExpr[[i]]$data <- multiExpr[[i]]$data[!outliers,]
  # multiMeta[[i]]$data <- multiMeta[[i]]$data[!outliers,]
}


## Network Construction

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
verbose = 5,networkType="signed")[[2]]);
# Plot the results:
pdf("../figures/softPower.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
for (col in 1:length(plotCols))
{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
}
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col]);
addGrid();
}
if (col==1)
{
text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
labels=powers,cex=cex1,col=colors[set]);
} else
text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set]);
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
} else
legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()



## cWGCNA Network

softPower=20
pairMat="AllStrains_CTX_6M"
cQuant=0.5 #consensusQuantile

net=blockwiseConsensusModules(multiExpr= multiExpr,maxBlockSize=25000,power=softPower,
                                  mergeCutHeight=0.25,nThreads=24,
                                  individualTOMFileNames = paste("../data_generated/signed_", pairMat,"TOM-Set%s-Block%b.rda",sep=""),
                                  corType="bicor",consensusQuantile=cQuant,
                                  minModuleSize=20,pamStage=FALSE,
                                  reassignThresholdPS = 1e-10,verbose=3,networkType="signed",
                                  saveTOMs = TRUE,
                                  consensusTOMFileNames = paste("../data_generated/signed_",
                                   pairMat,i,"consensusTOM-block.%b.rda",sep=""))

           load(file=paste("../data_generated/signed_", pairMat,i,"consensusTOM-block.1.rda",sep=""))
      geneTree= hclust(1-consTomDS,method="average");

   #### extracting goodGenes only
      multiExpr[[1]]$data=multiExpr[[1]]$data[,net$goodGenes]
 multiExpr[[2]]$data=multiExpr[[2]]$data[,net$goodGenes]
 multiExpr[[3]]$data=multiExpr[[3]]$data[,net$goodGenes]


 multiExpr=vector(mode="list",length=nSets)
multiExpr[[1]] = list(data=as.data.frame(datExpr_C57BL6))
names(multiExpr[[1]]$data)=colnames(datExpr_C57BL6)
rownames(multiExpr[[1]]$data)=rownames(datExpr_C57BL6)

multiExpr[[2]] = list(data=as.data.frame(datExpr_DBA))
names(multiExpr[[2]]$data)=colnames(datExpr_DBA)
rownames(multiExpr[[2]]$data)=rownames(datExpr_DBA)

multiExpr[[3]] = list(data=as.data.frame(datExpr_FVB))
names(multiExpr[[3]]$data)=colnames(datExpr_FVB)
rownames(multiExpr[[3]]$data)=rownames(datExpr_FVB)

checkSets(multiExpr) # check data size



       mColorh <- mLabelh <- colorLabels <- NULL
      for (minModSize in c(40,100,160)) {
        for (dthresh in c(0.1,0.2,0.25)) {
          for (ds in c(2,4)) {
            print("Trying parameters:")
            print(c(minModSize,dthresh,ds))
            tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
              minClusterSize = minModSize, cutHeight = 0.9999,
              deepSplit = ds, distM = as.matrix(1-consTomDS))

            merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                        cutHeight = dthresh)
            mColorh <- cbind(mColorh,labels2colors(merged$colors))
            mLabelh <- c(mLabelh,paste("DS=",ds," mms=",minModSize," dcor=",dthresh))
                }
        }
      }

      save(geneTree,mColorh,mLabelh,file=paste("./data/signed_",pairMat,i ,"_consensusTOMcuts.rda",sep=""))

       pdf(file=paste("./figures/signed_",pairMat, i ,"_consensus_TOMcuts.pdf",sep=""),height=8,width=11)
      plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters");
      dev.off()
      }


  rm(consTomDS)
save(list=ls(),file="consensus.rda")





# consMEs = net$multiMEs;
# moduleLabels = net$colors;
# # Convert the numeric labels to color labels
# moduleColors = labels2colors(moduleLabels)
# consTree = net$dendrograms[[1]];

# pdf("SignedDendro_Consensus.pdf",height=10, width=15)
# plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
# main = "Consensus gene dendrogram and module colors")
# dev.off()


targets_C57BL6=as.data.frame(targets[targets$Strain=="C57"|targets$Strain=="C57BL6",])
targets_DBA =as.data.frame(targets[targets$Strain=="DBA",])
targets_FVB =as.data.frame(targets[targets$Strain=="FVB",])

 multiMeta=vector(mode="list",length=nSets)

multiMeta[[1]] = list(data=as.data.frame(targets_C57BL6))
names(multiMeta[[1]]$data)=colnames(targets_C57BL6)
rownames(multiMeta[[1]]$data)=rownames(targets_C57BL6)

multiMeta[[2]] = list(data=as.data.frame(targets_DBA))
names(multiMeta[[2]]$data)=colnames(targets_DBA)
rownames(multiMeta[[2]]$data)=rownames(targets_DBA)

multiMeta[[3]] = list(data=as.data.frame(targets_FVB))
names(multiMeta[[3]]$data)=colnames(targets_FVB)
rownames(multiMeta[[3]]$data)=rownames(targets_FVB)

checkSets(multiMeta) # check data size

 tmpMulti =vector(mode="list",length=nSets)


####### Set 1 C57
thisMeta <- multiMeta[[1]]$data
  thisExpr <- t(multiExpr[[1]]$data)

  tmpMulti[[1]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"RIN"]))

  rownames(tmpMulti[[1]]$traitmat) <- rownames(thisMeta)
  colnames(tmpMulti[[1]]$traitmat) <- c("Wt.Tg","RIN")

  geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  ## Find adjusted multiple R^2 for each gene with each categorical variable
  for (i in 1:ncol(geneSigs)) {
    exprvec <- as.numeric(thisExpr[i,])
    # ager <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[1]]$traitmat[,"Age"])))$adj.r.squared,0))
    conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[1]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    rinr <- bicor(tmpMulti[[1]]$traitmat[,"RIN"],exprvec)
    geneSigs[,i] <- c(conditionr,rinr)
  }
  tmpMulti[[1]]$genesigs <- geneSigs

  # geneSigs[1,] <- numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  geneSigs[1,] =  numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
   geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  rownames(geneSigs) <- c("C57_Wt.Tg","C57_RIN")
  tmpMulti[[1]]$genecols <- geneSigs

  tmpMulti[[1]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  tmpMulti[[1]]$netData$TOMdendrogram <- geneTree
  tmpMulti[[1]]$netData$moduleColors <- mColorh
  tmpMulti[[1]]$netData$cutParameters <- mLabelh
  tmpMulti[[1]]$netData$annotColors <- geneSigs

  ## Second network - DBA
 thisMeta <- multiMeta[[2]]$data
  thisExpr <- t(multiExpr[[2]]$data)

  tmpMulti[[2]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"RIN"]))

  rownames(tmpMulti[[2]]$traitmat) <- rownames(thisMeta)
  colnames(tmpMulti[[2]]$traitmat) <- c("Wt.Tg","RIN")

  geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  ## Find adjusted multiple R^2 for each gene with each categorical variable
  for (i in 1:ncol(geneSigs)) {
    exprvec <- as.numeric(thisExpr[i,])
    # ager <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[2]]$traitmat[,"Age"])))$adj.r.squared,0))
    conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[2]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    rinr <- bicor(tmpMulti[[2]]$traitmat[,"RIN"],exprvec)
    geneSigs[,i] <- c(conditionr,rinr)
  }
  tmpMulti[[2]]$genesigs <- geneSigs

  # geneSigs[1,] <- numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  geneSigs[1,] =  numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
   geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  rownames(geneSigs) <- c("DBA_Wt.Tg","DBA_RIN")
  tmpMulti[[2]]$genecols <- geneSigs

  tmpMulti[[2]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  tmpMulti[[2]]$netData$TOMdendrogram <- geneTree
  tmpMulti[[2]]$netData$moduleColors <- mColorh
  tmpMulti[[2]]$netData$cutParameters <- mLabelh
  tmpMulti[[2]]$netData$annotColors <- geneSigs


    ############3rd network - FVB

 thisMeta <- multiMeta[[3]]$data
  thisExpr <- t(multiExpr[[3]]$data)

  tmpMulti[[3]]$traitmat <- cbind(as.factor(thisMeta[,"Wt.Tg"]),as.numeric(thisMeta[,"RIN"]))

  rownames(tmpMulti[[3]]$traitmat) <- rownames(thisMeta)
  colnames(tmpMulti[[3]]$traitmat) <- c("Wt.Tg","RIN")

  geneSigs <- matrix(NA,nrow=2,ncol=nrow(thisExpr))
  ## Find adjusted multiple R^2 for each gene with each categorical variable
  for (i in 1:ncol(geneSigs)) {
    exprvec <- as.numeric(thisExpr[i,])
    # ager <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[1]]$traitmat[,"Age"])))$adj.r.squared,0))
    conditionr <- sqrt(max(summary(lm(exprvec ~ as.factor(tmpMulti[[3]]$traitmat[,"Wt.Tg"])))$adj.r.squared,0))
    rinr <- bicor(tmpMulti[[3]]$traitmat[,"RIN"],exprvec)
    geneSigs[,i] <- c(conditionr,rinr)
  }
  tmpMulti[[3]]$genesigs <- geneSigs

  # geneSigs[1,] <- numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
  geneSigs[1,] =  numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
   geneSigs[2,] <- numbers2colors(as.numeric(geneSigs[2,]),blueWhiteRed(100),signed=TRUE,centered=TRUE,lim=c(-1,1))
  rownames(geneSigs) <- c("FVB_Wt.Tg","FVB_RIN")
  tmpMulti[[3]]$genecols <- geneSigs

  tmpMulti[[3]]$netData$netName <- c(paste("Signed bicor consensus network at quantile of 0.5 and power of 24"))
  tmpMulti[[3]]$netData$TOMdendrogram <- geneTree
  tmpMulti[[3]]$netData$moduleColors <- mColorh
  tmpMulti[[3]]$netData$cutParameters <- mLabelh
  tmpMulti[[3]]$netData$annotColors <- geneSigs



###################  plot

    mColorh1 <- cbind(mColorh,t(tmpMulti[[1]]$genecols),t(tmpMulti[[2]]$genecols),t(tmpMulti[[3]]$genecols))
  mLabelh1 <- c(mLabelh,rownames(tmpMulti[[1]]$genecols),rownames(tmpMulti[[2]]$genecols),rownames(tmpMulti[[3]]$genecols))
  i=0.5
       pdf(file=paste("./figures/signed_",pairMat, i ,"_consensus_TOMwithGenesSigs.pdf",sep=""),height=8,width=11)

  plotDendroAndColors(geneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Signed bicor consensus network at quantile of 0.5 and power of 24")
    dev.off()


  multiDataA_0.5 <- tmpMulti[[1]]
  multiDataB_0.5 <- tmpMulti[[2]]
  multiDataC_0.5 <- tmpMulti[[3]]
  save(multiDataA_0.5,multiDataB_0.5,multiDataC_0.5,file=paste("./data/",pairMat,"0.5_multiData.rda",sep=""))


		mms=100
		ds =4
		dthresh=0.2
        tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(1-consTomDS))

        merge <- mergeCloseModules(exprData = multiExpr,colors = tree$labels, cutHeight = dthresh)
        mColorh.cons <- labels2colors(merge$colors)
        mLabelh.cons <- "Merged Colors"
		mColorh1 <- cbind(mColorh.cons,t(multiDataA_0.5$genecols),t(multiDataB_0.5$genecols),t(multiDataC_0.5$genecols))
 		 mLabelh1 <- c(mLabelh.cons,rownames(multiDataA_0.5$genecols),rownames(multiDataB_0.5$genecols),rownames(multiDataC_0.5$genecols))
 		 i=0.5

     pdf(file=paste("./figures/signed_",pairMat, i ,"_consensus_TOMSelected_withGeneSigs.pdf",sep=""),height=8,width=11)
  plotDendroAndColors(geneTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh));
  dev.off()

mergedColors = labels2colors(merge$colors);
mergedMEs = merge$newMEs;
consensusMEs=rbind(mergedMEs[[1]]$data,mergedMEs[[2]]$data,mergedMEs[[3]]$data)
rownames(consensusMEs)=targets$Sample.ID

moduleColors = mergedColors

###change labels 2 colors
modCol=as.data.frame(table(mergedColors))
modLab=as.data.frame(table(merge$colors))
modLab$Var1=paste("ME",modLab$Var1,sep="")
mergeLab=merge(modLab, modCol,by="Freq")
indLab=match(rownames(moduleTraitPvalue),mergeLab$Var1)
rownames(moduleTraitPvalue)=as.character(mergeLab$mergedColors[indLab])
rownames(moduleTraitPvalue) -> rownames(moduleTraitCor)
rownames(moduleTraitPvalue) -> colnames(consensusMEs)

