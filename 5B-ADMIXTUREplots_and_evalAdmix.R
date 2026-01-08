library(reshape2)
library(ggplot2)

## READ IN DATA
setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/RADseq_Plate3_Results/admixture/')
sampIDs <- read.csv("../../SampleLocations/Marten_Hets_Fhom_coords.csv", fileEncoding="Latin1", header=TRUE, stringsAsFactors=FALSE)
prefix <- "11_LD50560" #high-missingness dataset
prefix <- "12_geno05mind10LD50560" #low-missingness dataset
prefix <- "11B_COMA_LD50560"#high-missingness with only COMA, ORcasc, and Lassen
prefix <- "11_LD50560noFirstRels"
prefix <- "11_LD_onesampperpop"
prefix <- "11B_COMAonly_LD50560"
prefix <- "11B_COMA_LD50560_noFirstRels"

sampinfo <- read.table(paste0(prefix,".fam"), sep=" ", header=FALSE)
k2 <- read.table(paste0(prefix, ".2.Q"))
k3 <- read.table(paste0(prefix, ".3.Q"))
k4 <- read.table(paste0(prefix, ".4.Q"))
k5 <- read.table(paste0(prefix, ".5.Q"))
k6 <- read.table(paste0(prefix, ".6.Q"))
k7 <- read.table(paste0(prefix, ".7.Q"))
k8 <- read.table(paste0(prefix, ".8.Q"))
k9 <- read.table(paste0(prefix, ".9.Q"))
k10 <- read.table(paste0(prefix, ".10.Q"))
k11 <- read.table(paste0(prefix, ".11.Q"))
k12 <- read.table(paste0(prefix, ".12.Q"))


## PLOT CV ERRORS AND LOG-LIKELIHOODS
#11_AB+LD50560
cv_a <- c(0.51032, 0.44301, 0.44044, 0.43240, 0.43136, 0.42708, 0.43307, 0.46383, 0.47432, 0.50033, 0.57600, 0.51032) 
logliks_a <- c( -544327.444039, -477663.048307, -454220.226845, -433718.617545, -419919.654071, -398106.481328, 
              -384268.702675, -372803.445570, -370398.560055, -359356.207708, -349527.503758, -345996.235939)
#11_AB+LD+relatives removed
cv_b <- c(0.55858,0.50353,0.51731,0.53125,0.55708,0.63384,0.70033,0.72694,.77630,0.89933,0.99782,0.99748) 
logliks_b <- c(-373066.477396,-324985.619836,-305062.93996,-292158.326410,-277565.623852,
               -267279.132705,-259340.369708,-250839.906425,-241774.483489,-235143.475009,-233206.393461,-225938.483981)

#11_AB+LD+onesampperpop
cv_c <- c(0.63219,0.68973,0.91204,1.09901,1.38043,1.53087,1.69740, 1.19547,1.95846,1.80406,1.43100,1.33853)
logliks_c <- c(-158371.225831,-140270.453540,-129989.921170,-123680.835955, -116998.123701,-110710.930613,-105678.700754,
               -100008.654370,-96958.883299,-94248.615557,-87014.633262,-85611.646007)


#12_geno5mind10_LD50560
cv <- c(0.48654,0.41629, 0.40737,0.39619, 0.39953, 0.39850, 0.44191, 0.45630, 0.44135, 0.47082, 0.54531, 0.58170)
logliks <- c( -132189.092360, -115256.381892, -110183.868040, -104110.110879, -101067.265308, -95607.132026, 
              -93135.006142, -91104.106668, -88377.532756, -86519.832828, -84313.337247, -83911.713398)

#11B COMA Only
cv <- c(0.51698, 0.48028, 0.46049, 0.46478, 0.46755, 0.54109, 0.57541, 0.60261, 0.73575)
logliks <- c(-352020.050266,-318263.853665,-297018.564904,-283554.440864,-270145.538130,-259631.703208,
             -252373.788378,-246184.696594,-241167.130869)


#11B_COMA/nearby only
cv <- c(0.51698,0.48028,0.46049,0.46478,0.46755,0.54109,0.57541,0.60261,0.73575,0.71712)
logliks <- c(-352020.050266,-318263.853665,-297018.564904,-283554.440864,-270145.538130,
             -259631.703208,-252373.788378,-246184.696594,-241167.130869, -235998.416169)
#12 LD with relatives removed
logliks <- c(-95859.092195,-83614.862399,-79257.913314,-74523.346296, -71594.131327,
             -69231.003737,-67412.014973,-65076.240627,-63995.320783,-61946.412671,
             -60989.462907,-60132.650012)


## CV error plot
par(mfrow=c(3,1), mar=c(2,2,1,1))
plot(1:12, cv_a, pch=16, xlab=NULL, ylab=NULL, main="All", xaxt="n", yaxt="n")
lines(1:12, cv_a, pch=16, lty="dashed")
points(2, cv_a[2], pch=8, cex=1.5)

plot(1:12, cv_b, pch=16, xlab=NULL, ylab=NULL,  main="No Close Kin", xaxt="n", yaxt="n")
lines(1:12, cv_b, pch=16, lty="dashed")
points(2, cv_b[2], pch=8, cex=1.5)

plot(1:12, cv_c, pch=16, xlab=NULL, ylab=NULL,  main="One Per Cluster", yaxt="n")
lines(1:12, cv_c, pch=16, lty="dashed")
points(2, cv_c[2], pch=8, cex=1.5)


## delta-logliks plots
dll_a <- logliks_a[2:12]-logliks_a[1:11]
dll_b <- logliks_b[2:12]-logliks_b[1:11]
dll_c <- logliks_c[2:12]-logliks_c[1:11]

plot(2:12, dll_a, pch=16, xlab=NULL, main="All", xaxt="n", yaxt="n")
lines(2:12, dll_a, lty="dashed")
points(2,dll_a[1], pch=8, cex=1.5)

plot(2:12, dll_b, pch=16, xlab=NULL, main="No Close Kin", xaxt="n", yaxt="n")
lines(2:12, dll_b, lty="dashed")
points(2,dll_b[1], pch=8, cex=1.5)

plot(2:12, dll_c, pch=16, xlab=NULL, main="One Per Cluster", yaxt="n")
lines(2:12, dll_c, lty="dashed")
points(2,dll_c[1], pch=8, cex=1.5)


# delta-logliks
plot(1:12, logliks, main="ADMIXTURE Loglikelihoods", xlab="K", ylab="Log-likelihood")

delta_logliks <- logliks[2:10]-logliks[1:9]
plot(2:10, delta_logliks, pch=16, xlab="K", ylab=expression(paste(Delta," log-likelihood")),
     main=expression(paste("ADMIXTURE ",Delta," log-likelihood")))
lines(2:10, delta_logliks, lty="dashed")
points(4, delta_logliks[3], pch=8, cex=1.5)


## SET INDIVIDUAL IDS
#inds <- sampinfo$V2
#inds <- unlist(lapply(inds, function(x){
#  split <- strsplit(x, "_")[[1]]
#  if(length(split)>1) return(paste(split[2],split[3],sep="_"))
#  if(length(split)==1) return(paste0(split[1],"LT"))
#}))#function
#inds <- paste(sampinfo$Pop, inds, sep="_")

inds <- sampIDs$MergeID[match(sampinfo$V2, paste("0",sampIDs$MergeIDraw,sampIDs$MergeIDraw, sep="_"))]
inds[which(inds=="M07B")] <- "SOreg_M07"
inds[c(4,27,32,44,58)] <- c("ORbmt_F10","SDunes_M6860","SOreg_M07","SOreg_F10","SDunes_Mrdkl")#12

inds <- factor(inds, levels=c("CO_M21256","CO_M21257","CO_F21258","CO_M21921","CO_M21922",
                              "WY_M31","WAcasc_M215","WAcasc_M40858",
                              "ORbmt_M01","ORbmt_03","ORbmt_M04","ORbmt_M06","ORbmt_M07",
                              "ORbmt_F08","ORbmt_M11","ORbmt_F10",
                              "Lssn_M61","Lssn_M64","Lssn_M65","Lssn_M66","Lssn_M67",
                              "ORcasc_F8943","ORcasc_F252","ORcasc_F258",
                              "Trinidad_M12","Trinidad_M01","Trinidad_F01","Trinidad_M21",
                              "NCali_Mrdkl1","NCali_Mrdkl2","NCali_Mrdkl3","NCali_M05","NCali_F03","NCali_M02","NCali_M03",
                              "SOreg_M15","SOreg_F05","SOreg_M13","SOreg_M08","SOreg_F10",
                              "SOreg_M07","SOreg_F11","SOreg_M17","SOreg_F04","SOreg_M16",
                              "SOreg_M09","SOreg_M10","SOreg_M12",
                              "SOreg_M22","SOreg_F07","SOreg_F08",
                              "SDunes_F04","SDunes_Mrdkl","SDunes_M6860","SDunes_M6856","SDunes_M04","SDunes_M01","SDunes_F05",
                              "NDunes_F07","NDunes_F01","NDunes_F02","NDunes_F06","NDunes_M02","NDunes_F31"))
 inds <- c("C21_215","D11_252","D11_258","BMT10","BMT11","BMT4","BMT6","BMT7","MACA01","Trinidad_F01",
          "NCali_F03","NCali_M03","NCali_M05","NDunes_F07","SDunes_M01","SDunes_M04", "WY_F31","AMMA_8943",
          "BMT8","SDunes_6856","SOreg_F05","SOreg_F07","SOreg_M09","SOreg_M13","NDunes_MACA31",
          "Lssn_M64","Lssn_M65","Lssn_M67","SOreg_F10","SOreg_F11","SOreg_M15","Trinidad_M21","SOreg_M22",
          "CraigsBeach","Newton2024","ZM21256","ZM21257","ZM21258","ZM21921","ZM21922","Charleston")
inds <- factor(inds, levels=c("ZM21256","ZM21257","ZM21258","ZM21921","ZM21922","WY_F31",
                              "BMT10","BMT11","BMT4","BMT6","BMT7","BMT8",
                              "C21_215","Lssn_M64","Lssn_M65","Lssn_M67","D11_252","D11_258","AMMA_8943",
                              "Trinidad_F01","Trinidad_M21","NCali_F03","NCali_M03","NCali_M05","CraigsBeach","Newton2024",
                              "SOreg_F05","SOreg_F07","SOreg_M09","SOreg_M13","SOreg_F10","SOreg_F11","SOreg_M15","SOreg_M22",
                              "SDunes_M01","SDunes_M04","SDunes_6856","Charleston",
                              "MACA01","NDunes_F07","NDunes_MACA31"))


k2$inds <- k3$inds <- k4$inds <- k5$inds <- k6$inds <- k7$inds <- k8$inds <- k9$inds <- inds



### PLOTS!
orig <- par('mar')
par(mar=c(10,4,4,1), mfrow=c(4,1))
cols <- rainbow(6)


## K2 PLOT
#k2e <- t(as.matrix(k2))
#k2e <- k2e[,order(k2e[1,])]
#barplot(k2e[-which(rownames(k2e)=="inds"),], names.arg=k2e[which(rownames(k2e)=="inds"),], 
#        ylab="Ancestry", col=cols[c(1,3)], las=2, main=paste("K=2,",expression("delta")),
#        cex.names=0.9, cex.axis=0.7, space=rep(0,ncol(k2e)), border=NA)

k2long <- reshape2::melt(k2, id.vars="inds", value.name="Q")
k2long$inds <- factor(k2long$inds, levels=levels(inds))
k2p <- ggplot(k2long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=cols[c(5,1)])+
  ylab(element_blank())+
  xlab(element_blank())


## K3 PLOT
k3long <- reshape2::melt(k3, id.vars="inds", value.name="Q")
k3long$inds <- factor(k3long$inds, levels=levels(inds))
k3p <- ggplot(k3long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c(cols[1],cols[6],cols[5]))+#11E 
  ylab(element_blank())+
  xlab(element_blank())


## K4 PLOT
k4long <- reshape2::melt(k4, id.vars="inds", value.name="Q")
k4long$inds <- factor(k4long$inds, levels=levels(inds))
k4p <- ggplot(k4long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=cols[c(4,6,1,5)])+#11LD one per pop
  #scale_fill_manual(values=cols[c(1,2,4,5)])+#11LD no rels
  #scale_fill_manual(values=cols[c(1,4,5,2)])+#12+#1,2,4,5 #11AB_LD
  #scale_fill_manual(values=c(cols[2],cols[5],"purple",cols[c(1)]))+#11B_COMA
  ylab(element_blank())+
  xlab(element_blank())


## K5 PLOT
k5long <- reshape2::melt(k5, id.vars="inds", value.name="Q")
k5long$inds <- factor(k5long$inds, levels=levels(inds))
k5p <- ggplot(k5long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  #scale_fill_manual(values=c(cols[c(2,5,4,1)],"orange"))+
  #scale_fill_manual(values=c(cols[c(1,4,5,2)],"orange"))+#cols[c(4,5)],"orange",cols[c(1,2)])#11e
  #scale_fill_manual(values=c(cols[2],cols[1],cols[5],cols[4],"purple"))+#11B_COMA
  #scale_fill_manual(values=c(cols[4],cols[6],cols[1],cols[5],cols[2]))+#11LD_noRels
  scale_fill_manual(values=c("orange",cols[4],cols[5],cols[1],cols[6]))+#11LD_oneSamp
  #scale_fill_manual(values=c(cols[2],cols[1],"orange",cols[5],"purple"))+#11B_COMAonly
  ylab(element_blank())+
  xlab(element_blank())+
  theme(plot.margin=unit(c(0,0,0.2,1),'lines'))


## K6 PLOT
k6long <- reshape2::melt(k6, id.vars="inds", value.name="Q")
k6long$inds <- factor(k6long$inds, levels=levels(inds))
k6p <- ggplot(k6long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  #scale_fill_manual(values=c(cols[c(2,1,6,4,5)],"purple"))+#11
  #scale_fill_manual(values=c(cols[c(6,5,2,1)],"orange",cols[4]))+#12
  #scale_fill_manual(values=c(cols[6],cols[1],"orange",cols[2],cols[5],cols[4]))+#11LD_noRels
  scale_fill_manual(values=c(cols[2],"orange",cols[5],cols[6],cols[1],"sienna"))+#11LD_onesampperpop
  #scale_fill_manual(values=c(cols[1],"orange",cols[4],"purple",cols[2],cols[5]))+#11B_COMAonly
  ylab(element_blank())+
  xlab(element_blank())


## K7 PLOT
k7long <- reshape2::melt(k7, id.vars="inds", value.name="Q")
k7long$inds <- factor(k7long$inds, levels=levels(inds))
#k7p <- 
ggplot(k7long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  #scale_fill_manual(values=c("purple",cols[c(6,1,5)],cols[2],cols[4],cols[3]))+#11
  #scale_fill_manual(values=c(cols[c(6,5,2,1)],"orange",cols[4]))+#12
  #scale_fill_manual(values=c("orange",cols[4],cols[6],"purple",cols[3],cols[2],cols[1]))+#11LD_noRels
  scale_fill_manual(values=c(cols[6],cols[1],cols[2],"purple",cols[4],"sienna",cols[5]))+
  ylab(element_blank())+
  xlab(element_blank())

k8long <- reshape2::melt(k8, id.vars="inds", value.name="Q")
k8long$inds <- factor(k8long$inds, levels=levels(inds))
k8p <- ggplot(k8long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  #scale_fill_manual(values=c(cols[5],"orange",cols[6],"purple",cols[1],cols[4],cols[2],cols[3]))+#11
  #scale_fill_manual(values=c(cols[c(6,5,2,1)],"orange",cols[4]))+#12
  scale_fill_manual(values=c(cols[1],"orange",cols[3],cols[4],cols[2],cols[6],"purple",cols[5]))+#11LD_noRels
  ylab(element_blank())+
  xlab(element_blank())


k9long <- reshape2::melt(k9, id.vars="inds", value.name="Q")
k9long$inds <- factor(k9long$inds, levels=levels(inds))
k9p <- ggplot(k9long)+
  geom_bar(aes(x=inds, y=Q, fill=variable), stat="identity", width=1)+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  #scale_fill_manual(values=c(cols[2],cols[4],"turquoise",cols[3],cols[1],cols[6],"orange","purple",cols[5]))+#11
  #scale_fill_manual(values=c(cols[c(6,5,2,1)],"orange",cols[4]))+#12
  scale_fill_manual(values=c("sienna","orange",cols[3], cols[2],cols[5],"purple",cols[4],cols[1], cols[6]))+#11LD_noRels
  ylab(element_blank())+
  xlab(element_blank())+
  theme(plot.margin=unit(c(0,0,0.1,0),'lines'))


k2p <- k2p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k3p <- k3p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k4p <- k4p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k5p <- k5p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k6p <- k6p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k7p <- k7p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k8p <- k8p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
k9p <- k9p+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(plot.margin=unit(c(0,0,0.1,0),'lines'))
ggarrange(k2p,k3p,k4p,k5p,
          #,k6p,k7p,k8p,k9p,
          nrow=4,ncol=1)


p5 <- read.table(paste0(prefix,".5.P"))
ggplot(p5)+
  geom_violin(aes(x=1,y=V1), col=cols[4])+
  geom_violin(aes(x=2,y=V2), col="orange")+
  geom_violin(aes(x=3,y=V3), col=cols[5])+
  geom_violin(aes(x=4,y=V4), col=cols[2])+
  geom_violin(aes(x=5,y=V5), col=cols[1])
ggplot(p6)+
  geom_violin(aes(x=1,y=V1), col=cols[4])+
  geom_violin(aes(x=2,y=V2), col="orange")+
  geom_violin(aes(x=3,y=V3), col=cols[5])+
  geom_violin(aes(x=4,y=V4), col=cols[2])+
  geom_violin(aes(x=5,y=V5), col=cols[1])+
  geom_violin(aes(x=6,y=V6), col="purple")

###########################
#### evalAdmix outputs ####
###########################
source("evalAdmix/visFuns.R")

# read in individuals and set populations
prefix <- "11_LD50560noFirstRels"
sampinfo <- read.table(paste0(prefix,".fam"))
inds <- sampIDs$MergeID[match(sampinfo$V2, paste("0", sampIDs$MergeIDraw,sampIDs$MergeIDraw, sep="_"))]
pop <- unlist(lapply(inds, function(x) strsplit(x,"_")[[1]][1]))
#pop <- factor(pop, levels=c("CO","WY","ORbmt","WAcasc","Lssn","ORcasc","Trinidad","NCali","SOreg","SDunes","NDunes"))

# read in eval file and admixture props
k <- 9
evalk2 <- read.table(paste0(prefix,".evalADMIXK",k,".out"))
q <- read.table(paste0(prefix,".",k,".Q"))

# plot, order by pop
ord <- orderInds(pop=pop, q=q, popord=c("CO","WY","WAcasc","ORbmt","Lssn","ORcasc","NCali","Trinidad","SOreg","SDunes","NDunes"))
plotCorRes(cor_mat=evalk2, ord=ord, pop=pop, max_z=0.4, min_z=-0.4, 
           #title="K=2", #title=paste0("Evaluation of admixture proportions with K=",as.character(k)), 
           rotatelabpop=90, rotatelabsuperpop=1, pop_labels=c(F,T), adjlab=0.2)


## aggregate for each group for mapping
k4$pop <- unlist(lapply(inds, FUN=function(x) strsplit(x, "_")[[1]][1]))
unique(k4$pop)
qbypop <- aggregate(data=k4, V1~pop, FUN=sum)
qbypop$V2 <- aggregate(data=k4, V2~pop, FUN=sum)[,2]
qbypop$V3 <- aggregate(data=k4, V3~pop, FUN=sum)[,2]
qbypop$V4 <- aggregate(data=k4, V4~pop, FUN=sum)[,2]
head(qbypop)
write.csv(qbypop, paste0(prefix,"_ByPop.Q.csv"))
