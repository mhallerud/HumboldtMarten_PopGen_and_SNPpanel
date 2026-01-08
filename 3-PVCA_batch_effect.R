#-------------PRINCIPAL VARIANCE COMPONENTS ANALYSIS FOR BATCH EFFECTS-----------------#
# Author: Maggie Hallerud
# Paper: Humboldt marten ddRADseq
setwd("/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/RADseq_Plate3_Results/VCFs/")

# load dependencies
#library(pvca)
library(lme4)
library(vcfR)
library(adegenet)
library(ggplot2)



########################################
#### RUN PVCA FOR FILTERED DATASET #####
########################################
#### RUN PVCA FOR MINIMALLY FILTERED DATASET ####
# read in minimally filtered data
vcf <- read.vcfR("4_populations.snps.Q20.vcf")
genind <- vcfR2genind(vcf)
rawtab <- tab(genind, freq=T, NA.method="mean")

rawbatch <- c(rep("Plate1",40),rep("Plate2",26),rep("Plate3",29))
cbind(indNames(genind), rawbatch)#check!

# extract eigenvalues and compute eigenvectors from loadings
eigens <- eigen(cor(t(rawtab)))
rawvals <- eigens$values
rawvecs <- eigens$vectors

# set up data frame to hold PVCA output- rows=PCs, cols=variation per group
sum(rawvals[1:10])/sum(rawvals)#use the first 10 PCs- 93% of variation...
rawREM <- as.data.frame(matrix(NA, nrow=10, ncol=4))

names(rawREM) <- c("Batch1","Batch2","Batch3","Sigma")
batch1 <- batch2 <- batch3 <- rawbatch
batch1[batch1=="Plate1"] <- 1; batch1[batch1!=1] <- 0
batch2[batch2=="Plate2"] <- 1; batch2[batch2!=1] <- 0
batch3[batch3=="Plate3"] <- 1; batch3[batch3!=1] <- 0
rownames(rawREM) <- paste0("PC",1:10)

# run mixed-effects model with batch as random intercept and export associated variances for PCs 1-10
for (pc in 1:10){
  # make dataframe with PC values in first column and batch dummy var in next 3 cols
  df <- data.frame(PC=rawvecs[,pc], Batch1=as.factor(batch1), Batch2=as.factor(batch2), Batch3=as.factor(batch3))
  # run LMER on PC with batches as dummy variables
  reml <- lme4::lmer(PC ~ (1|Batch1)+(1|Batch2)+(1|Batch3), data=df, REML=TRUE, verbose=FALSE, na.action=na.omit)
  # extract variance and correlation components of model
  rawREM$Batch1[pc] <- unlist(VarCorr(reml))[1]#variance explained for batch1
  rawREM$Batch2[pc] <- unlist(VarCorr(reml))[2]#variance explained for batch2
  rawREM$Batch3[pc] <- unlist(VarCorr(reml))[3]#variance explained for batch3
  rawREM$Sigma[pc] <- sigma(reml)^2 #residual standard deviation
  # convert to % model variance explained per batch
  rawREM[pc,] <- rawREM[pc,] / sum(rawREM[pc,])
  # weight by % covariance explained based on eigenvals
  #rawREM[pc,] <- rawREM[pc,] * (rawvals[pc]/sum(rawvals))
}#pc

# plot contributions per PC and batch effects per PC
rawREM$EigenVal <- rawvals[1:10]/sum(rawvals)
rawREM$PC <- factor(rownames(rawREM), levels=paste0("PC",1:10))
rawREM_long <- reshape(rawREM, varying=names(rawREM)[1:5], times=c("Batch1","Batch2","Batch3","Sigma","EigenVal"), 
                       idvar=c("PC"), timevar="Effect", v.names="Estimate", direction="long")
rawREM_long$Effect[rawREM_long$Effect=="Sigma"] <- "Residual"
ggplot()+
  #geom_bar(data=rawREM_long[rawREM_long$Effect!="Sigma",], aes(x=PC, y=Estimate), stat="identity", fill="gray80")+
  #geom_bar(data=rawREM_long[rawREM_long$Effect=="EigenVal",], aes(x=PC, y=Estimate), stat="identity", fill="gray20")+
  geom_bar(data=rawREM_long[rawREM_long$Effect!="EigenVal",], aes(x=PC, y=Estimate, fill=Effect), stat="identity")+
  theme_classic()+
  labs(fill="Effect")+
  xlab("Principal Component")+
  ylab("Model variance")+
  scale_y_continuous(labels=scales::percent, expand=c(0,0))+
  scale_fill_manual(values=c("red","gold2","blue","gray20"))+
  ggtitle("Batch Effects per PC Axis")+
  theme(text=element_text(size=12), axis.text.x=element_text(angle=90), title=element_text(face="bold"), legend.position="none")

# convert to % total variation of PCA and plot
rawREM[,1:4]/rawREM$EigenVal #scale by eigenval
rawREMsum <- colSums(rawREM[,1:4]) / sum(rawREM[,1:4]); rawREMsum
names(rawREMsum)[4] <- "Residual"
barplot(rawREMsum,  xlab=NULL,
        ylab = "Cumulative variance in PCA", ylim= c(0,1.1),
        col = c("red","gold2","blue","gray20"), las=2, main="Variation Explained in Raw Data")



#### RUN PVCA FOR FILTERED DATASET ####
# read in high-missingness dataset (includes replicates)
vcf <- read.vcfR("HighMissingnessRAD_withDups.vcf")
genind <- vcfR2genind(vcf)
filttab <- tab(genind, freq=F, NA.method="mean")

filtbatch <- c(rep("Plate1",29),rep("Plate2",22),rep("Plate3",24))
cbind(indNames(genind), filtbatch)#check!

# extract eigenvalues and compute eigenvectors from loadings
#filtpca <- dudi.pca(tab, scannf=FALSE, nf=25)
eigens <- eigen(cor(t(filttab)))
filtvals <- eigens$values
filtvecs <- eigens$vectors

# run mixed-effects model with batch as random intercept and export associated variances
sum(filtvals[1:10])/sum(filtvals)#use the first 10 PCs...
filtREM <- as.data.frame(matrix(NA, nrow=10, ncol=4))
names(filtREM) <- c("Batch1","Batch2","Batch3","Sigma")
batch1 <- batch2 <- batch3 <- filtbatch
batch1[batch1=="Plate1"] <- 1; batch1[batch1!=1] <- 0
batch2[batch2=="Plate2"] <- 1; batch2[batch2!=1] <- 0
batch3[batch3=="Plate3"] <- 1; batch3[batch3!=1] <- 0
for (pc in 1:10){
  df <- data.frame(PC=filtvecs[,pc], Batch1=as.factor(batch1), Batch2=as.factor(batch2), Batch3=as.factor(batch3))
  reml <- lmer(PC ~ (1|Batch1)+(1|Batch2)+(1|Batch3), data=df, REML=TRUE, verbose=FALSE, na.action=na.omit)
  filtREM$Batch1[pc] <- unlist(VarCorr(reml))[1]
  filtREM$Batch2[pc] <- unlist(VarCorr(reml))[2]
  filtREM$Batch3[pc] <- unlist(VarCorr(reml))[3]
  filtREM$Sigma[pc] <- sigma(reml)^2
  # convert to % model variance explained per PC
  filtREM[pc,] <- filtREM[pc,] / sum(filtREM[pc,])
  # weight by % covariance explained based on eigenvals
  #filtREM[pc,] <- filtREM[pc,] * (filtvals[pc]/sum(filtvals))
}#pc

# plot contributions per PC and batch effects per PC
filtREM$EigenVal <- filtvals[1:10]/sum(filtvals)
filtREM$PC <- factor(paste0("PC",1:10), levels=paste0("PC",1:10))
filtREM_long <- reshape(filtREM, varying=names(filtREM)[1:5], times=c("Batch1","Batch2","Batch3","Sigma","EigenVal"), 
                       idvar=c("PC"), timevar="Effect", v.names="Estimate", direction="long")
filtREM_long$Effect[filtREM_long$Effect=="Sigma"] <- "Residual"
ggplot()+
  #geom_bar(data=rawREM_long[rawREM_long$Effect!="Sigma",], aes(x=PC, y=Estimate), stat="identity", fill="gray80")+
  #geom_bar(data=rawREM_long[rawREM_long$Effect=="EigenVal",], aes(x=PC, y=Estimate), stat="identity", fill="gray20")+
  geom_bar(data=filtREM_long[filtREM_long$Effect!="EigenVal",], aes(x=PC, y=Estimate, fill=Effect), stat="identity")+
  theme_classic()+
  labs(fill="Effect")+
  xlab("Principal Component")+
  ylab("Model variance")+
  scale_y_continuous(labels=scales::percent, expand=c(0,0))+
  scale_fill_manual(values=c("red","gold2","blue","gray20"))+
  ggtitle("Batch Effects per PC Axis")+
  theme(text=element_text(size=12), axis.text.x=element_text(angle=90), legend.position="none", title=element_text(face="bold"))

# convert to % total variation & plot
filtREM <- filtREM[,1:4]/filtREM$EigenVal #weight by eigenval
filtREMsum <- colSums(filtREM[,1:4]) / sum(filtREM[,1:4])
sum(filtREMsum[1:3]) #13% variation explained by batch effect
names(filtREMsum)[4] <- "Residual"
barplot(filtREMsum,  xlab=NULL,
        ylab = "Cumulative variance in PCA", ylim= c(0,1.1),
        col = c("red","gold2","blue","gray20"), las=2, main="Variation Explained in Filtered Data")
