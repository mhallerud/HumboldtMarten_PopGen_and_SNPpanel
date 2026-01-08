#-------------PRINCIPAL COORDINATES ANALYSIS FOR ddRADseq-----------------#
# Author: Maggie Hallerud
# Paper: Humboldt marten ddRADseq
setwd("/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/Marten_Population_Genetics_Code/data/")

# load dependencies
library(vcfR)
library(adegenet)
library(ggplot2)
library(ggpubr)
library(viridis)
library(scales)

# Load sample data
popmap <- read.table("COMBINED_marten_popmap.txt",col.names=c("ID","pop"))
popmap$pop[popmap$pop%in%c("Border","Redwoods","Trinidad")] <- "NCali"




######################
#### RADSEQ PCOA  ####
######################
# load VCF without technical replicates
finalvcf <- read.vcfR("HighMissingnessRAD.vcf")
finalgenind <- vcfR2genind(finalvcf)

# add pops based on IDs
pops <- pop(finalgenind) <- unlist(lapply(indNames(finalgenind), function(x) strsplit(x,"_")[[1]][1]))

# create population labels
pops[pops=="NDunes"] <- "N Dunes"
pops[pops=="SDunes"] <- "S Dunes"
pops[pops=="MapleCrk"] <- "Maple Creek"
pops[pops=="NCali"] <- "N California"
pops[pops=="WAcasc"] <- "WA Cascades"
pops[pops=="ORcasc"] <- "OR Cascades"
pops[pops=="Lssn"] <- "Lassen"
pops[pops=="ORbmt"] <- "Blue Mtns"
pops[pops=="WY"] <- "Wyoming"
pops[pops=="SOreg"] <- "S Oregon"
pops[pops=="CO"] <- "Colorado"

# factor (include unrepresented pops so that colors match GBAS results)
pops <- factor(pops, levels=c("N Coast","N Dunes","S Dunes","S Oregon","S Coast",
                              "Border","Maple Creek","N California","Applegate",
                              "OR Cascades","Lassen","WA Cascades","Blue Mtns",
                              "Wyoming","Colorado","Unk"))

# plate information
finalplate <- c(rep("Plate1",26), rep("Plate2",16), rep("Plate3",21))#10AB

# run principal coordinates analysis (PCO)
finaltab <- tab(finalgenind, freq=T, NA.method="mean")
finalpco <- dudi.pco(dist(finaltab), nf=50, scannf=FALSE)
summary(finalpco)

# default plots
s.class(finalpco$li, fac=pops, cstar=0, axesell=FALSE, col=rainbow(16))
s.class(finalpco$li, fac=as.factor(finalplate), cstar=0, axesell=FALSE, col=rainbow(3))
add.scatter.eig(finalpco$eig, posi="topright", xax=1, yax=2)

# color palette
pal <- c("#000000",#black
         "#db6d00",#orange
         "gray20",#"#252525",#darkgray
         "gray60",#676767",#lightgray
         "gray90",#"#ffffff",#white
         "mediumpurple4",#171723",#navy
         "darkcyan",#"#004949",#dark teal
         "aquamarine1",#009999",#light teal
         "#22cf22",#bright green
         "darkviolet",#"#490092",#purple
         "steelblue1",#"#006ddb",#sky blue
         "#b66dff",#light purple
         "#ff6db6",#bright pink
         "#920000",#dark red
         "#ffdf4d",#yellow
         "#8f4e00")#brown

# final plot (FIGURE 1 IN MS)
#(rotate axis 1 to align with geography)
ggplot(finalpco$li, aes(x=-A1, y=A2, col=pops))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  geom_point(aes(pch=finalplate), size=6, alpha=0.8)+
  #stat_ellipse(level=0.90)+
  theme_classic()+
  xlab("PCO1 (17.0%)")+
  ylab("PCO2 (8.0%)")+
  #labs(col="Population",pch="Plate")+
  scale_color_manual(breaks=levels(pops), values=rev(pal))+#pals::watlington(n=16)+#rainbow(16))+
  theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.line=element_blank(), axis.title=element_text(size=20), #panel.background = element_rect(fill="gray90"),
        plot.background=element_rect(fill=NA, color='black'), plot.title=element_text(face="bold",size=20))+
  coord_fixed()+
  ggtitle("ddRADseq (N=63, SNPs=14,652)")

## plot by missingness
finaltab <- tab(finalgenind, freq=T, NA.method="asis")
finalmiss <- apply(finaltab, 1, function(x) sum(is.na(x))/length(x))
missplot <- ggplot(finalpco$li, aes(x=-A1, y=A2, col=finalmiss))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  theme_classic()+
  geom_point(aes(pch=finalplate), size=3, alpha=0.7)+
  xlab("PCO1 (17.0%)")+
  ylab("PCO2 (8.0%)")+
  labs(col="Missingness",pch="Plate")+
  theme(axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.line=element_blank())+
  scale_color_viridis(option="plasma")+
  ggtitle("A")+  
  theme(panel.border = element_rect(color="black",fill=NA)) +
  coord_fixed()


## plot by sequencing depth
# depth file created by running "vcftools --vcf HighMissingnessRAD.vcf --depth"
depths <- read.table("out.idepth", header=T)
# check...
indNames(finalgenind)[!indNames(finalgenind)%in%depths$INDV]
sum(depths$INDV!=indNames(finalgenind))
# plot
depthplot <- ggplot(finalpco$li, aes(x=-A1, y=A2, col=depths$MEAN_DEPTH))+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  theme_classic()+
  geom_point(aes(pch=finalplate), size=3, alpha=0.7)+
  xlab("PCO1 (17.0%)")+
  ylab("PCO2 (8.0%)")+
  labs(col="Mean Depth",pch="Plate")+
  theme(axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.line=element_blank())+
  scale_color_viridis()+
  ggtitle("B")+
  coord_fixed()+
  theme(panel.border = element_rect(color="black",fill=NA))
  

## screeplot
screeplot(finalpco, main=NULL)
scree <- ggplot()+
  geom_bar(aes(x=1:length(finalpco$eig), y=finalpco$eig), stat="identity")+
  xlab("Axis")+ylab("Inertia")+#ggtitle("Screeplot")+
  theme_classic()+
  theme(text=element_text(size=12))

variance <- finalpco$eig / sum(finalpco$eig)
sumvariance <- unlist(lapply(1:length(variance), function(x) sum(variance[1:x])))
scree <- ggplot(mapping=aes(x=1:length(finalpco$eig), y=variance*100))+
  geom_bar(stat="identity", width=1)+
  #geom_line(lty="dashed")+
  theme_classic()+
  xlab("Axis")+ylab("% Variation explained")+
  theme(text=element_text(size=12))+
  ggtitle("D")+
  theme(panel.border = element_rect(color="black",fill=NA))


## plot other axes
s.class(finalpco$li, fac=pops, xax=2, yax=3, cstar=0, axesell=F, col=rainbow(16))
s.class(finalpco$li, fac=pops, xax=3, yax=4, cstar=0, axesell=F, col=rainbow(16))
pco23 <- ggplot(finalpco$li)+
  geom_point(aes(x=A2, y=A3, col=pops))+
  theme_bw()+
  scale_color_manual(values=rainbow(10))+
  xlab("PCO2 (8.0%)")+ylab("PCO3 (6.5%)")+
  ggtitle("C")+
  theme(legend.position="none")
pco34 <- ggplot(finalpco$li)+
  geom_point(aes(x=A3, y=A4), size=2.1)+#outline for pts
  geom_point(aes(x=A3, y=A4, col=pops),size=2)+
  theme_bw()+
  scale_color_manual(breaks=levels(pops), values=rev(pal))+#pals::watlington(n=16)+#rainbow(16))+
  xlab("PCO3 (6.5%)")+ylab("PCO4 (5.3%)")+
  ggtitle("C")+
  labs(col="Population")+
  coord_fixed()+
  theme(panel.border = element_rect(color="black",fill=NA), legend.title=element_blank())


ggpubr::ggarrange(missplot, depthplot, pco34, scree)



############################
#### HUMBOLDT-ONLY PCOA ####
############################
# subset data
comasub <- finalgenind[which(pop(finalgenind)%in%c("NDunes","SDunes","SOreg","Border","NCali","MapleCrk")),]
comatab <- tab(comasub)

# remove monomorphic sites- these won't affect PCO but to get count of SNPs
library(dartR)
comagl <- dartR::gl.filter.monomorphs(dartR::gi2gl(comasub))
comagl

# run PCO
comapco <- dudi.pco(dist(comatab), nf=50, scannf=FALSE)
summary(comapco)
#s.class(comapco$li, fac=pop(comasub), axesell=FALSE, cstar=0, col=rainbow(16))

# fix pops for plotting colors...
comapops <- as.character(pop(comasub))
comapops[comapops=="NDunes"] <- "N Dunes"
comapops[comapops=="SDunes"] <- "S Dunes"
comapops[comapops=="SOreg"] <- "S Oregon"
comapops[comapops=="NCali"] <- "N California"
comapops[comapops=="MapleCrk"] <- "Maple Creek"
cols <- rev(pal)[which(levels(pops) %in% comapops)]
comapops <- factor(comapops, levels=c("N Dunes","S Dunes","S Oregon","Border","Maple Creek","N California"))

# plot
ggplot()+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  geom_point(aes(x=comapco$li[,1], y=-comapco$li[,2], col=comapops), size=6, alpha=0.8)+
  xlab("PCO1 (15.3%)")+
  ylab("PCO2 (9.0%)")+
  #stat_ellipse(level=0.90)+
  theme_classic()+
  scale_color_manual(breaks=levels(comapops), values=cols)+
  theme(axis.text=element_blank(), axis.ticks=element_blank(), legend.position="none",  
        axis.line=element_blank(), text=element_text(size=16), #panel.background = element_rect(fill="gray90"),
        plot.background=element_rect(fill=NA, color='black'), plot.title=element_text(face="bold"))+
  coord_fixed()+
  ggtitle("ddRADseq (N=39, SNPs=10,206)")



####################################
#### MAPPING PCO AXES VS. SPACE ####
####################################
## read in coordinates
popmap <- read.csv("../SampleLocations/Sample_coords.csv", stringsAsFactors=FALSE, fileEncoding='Latin1', header=T)

# check for missing names
indNames(finalgenind)[!indNames(finalgenind)%in%popmap$MergeID]

# set up coordinate df
gtcoords <- data.frame(Longitude=popmap$Longitude[match(indNames(finalgenind), popmap$MergeID)],
                       Latitude=popmap$Latitude[match(indNames(finalgenind), popmap$MergeID)])
gtcoords$Longitude <- as.numeric(gtcoords$Longitude)
sum(is.na(gtcoords$Latitude))

# convert lat/long coords to a projection that preserves distance
library(sp)
crs <- CRS("EPSG:4326")
sp <- sp::SpatialPointsDataFrame(gtcoords[,c("Longitude","Latitude")], gtcoords, proj4string=crs)
sp <- spTransform(sp, CRS("ESRI:102005"))#USA contiguous equidistant conic

# extract transformed coordinates
gtcoords <- as.data.frame(coordinates(sp))
names(gtcoords) <- c("x","y")

# plot
ggplot()+
  geom_point(aes(x=coords$x, y=coords$y, col=pco$li[,1]), alpha=0.5, size=3)+
  theme_bw()+
  labs(col="PCO1")+
  ggtitle("PCo1 in Geographic Space")+
  coord_fixed(ratio=1)+
  theme(text=element_text(size=16), legend.position="none", axis.text=element_blank(), axis.title=element_blank())+
  scale_color_viridis()



###########################################
#### PCoAs AT VARIOUS FILTERING STAGES ####
###########################################
rawvcf <- read.vcfR("VCFs/populations.snps.vcf") #97 x 516545
vcf4 <- read.vcfR("VCFs/4_populations.snps.Q20.vcf") #97 x 75749
vcf5a <- read.vcfR("VCFs/5A_removeOutliers.recode.vcf") #94 x 46,433
vcf5b <- read.vcfR("VCFs/5B_geno30.vcf") #94 x 15762
vcf5c <- read.vcfR("VCFs/5C_mind50.vcf") #75 x 15299
vcf10 <- read.vcfR("VCFs/HighMissingnessRAD.vcf") #63 x 14,652 - duplicates removed!
vcf11 <- read.vcfR("VCFs/LDprunedRAD.vcf") #63 x 12,389
vcf12 <- read.vcfR("VCFs/LowMissingnessRAD.vcf") #58x3586
vcf12b <- read.vcfR("VCFs/12_geno05mind10LD50560.recode.vcf") #58 x 3055


#### 0: Raw Dataset ####
# extract tab
rawgenind <- vcfR2genind(rawvcf)
colnames(rawvcf@gt)[!colnames(rawvcf@gt)%in%indNames(rawgenind)]#MACA_HU_F05A
rawtab <- tab(rawgenind, freq=T, NA.method="asis")

# grab ids, plates, missingness
pop(rawgenind) <- popmap$pop[match(indNames(rawgenind), popmap$ID)]
pop(rawgenind) <- factor(pop(rawgenind), levels=c("WaCascades","OrCascades","Lassen","BlueMountains",
                                                  "SDunes","NDunes","NCali","Rockies","SouthernOregon",
                                                  "Unk","Denver"))
rawplate <- c(rep("Plate1",40),rep("Plate2",26),rep("Plate3",30))
rawmiss <- apply(rawtab, 1, function(x) sum(is.na(x))/length(x))

# run pco
rawpco <- dudi.pco(dist(rawtab), scannf=FALSE, nf=50)
summary(rawpco)

# plots
rawsub <- which(rawpco$li$A2< -50 | rawpco$li$A1>100)
indNames(rawgenind)[rawsub]
rawplot <- ggplot(rawpco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (9.5%)")+
  ylab("PCO2 (7.1%)")+
  ggtitle("A")
p1 <- rawplot+
  #ylab("PCO2 (7.1%)")+
  geom_point(aes(col=rawplate))+
  labs(color="Plate")
p3 <- rawplot+
  geom_point(aes(col=pop(rawgenind)))+
  labs(color="Population")
p2 <- rawplot+
  #xlab("PCO1 (9.5%)")+
  geom_point(aes(col=rawmiss))+
  geom_text(data=rawpco$li[rawsub,], aes(label=c("Lssn_M11","BMT_M05","NDunes_F03","","")))+
  labs(color="Missingness")+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p1,p2,p3, nrow=1)



#### 4: Minimum quality-filtered dataset ####
# extract tab
qualgenind <- vcfR2genind(vcf4)
colnames(vcf4@gt)[!colnames(vcf4@gt)%in%indNames(qualgenind)]#MACA_HU_F05A, MACAtail
qualtab <- tab(qualgenind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(qualgenind) <- popmap$pop[match(indNames(qualgenind), popmap$ID)]
qualplate <- c(rep("Plate1",40),rep("Plate2",26),rep("Plate3",29))
qualmiss <- apply(qualtab, 1, function(x) sum(is.na(x))/length(x))

# run pco
qualpco <- dudi.pco(dist(qualtab), scannf=FALSE, nf=50)
summary(qualpco)

# plots
qualsub <- which(qualpco$li$A2< -50)
indNames(qualgenind)[qualsub]
qualplot <- ggplot(qualpco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (10.7%)")+
  ylab("PCO2 (12.3%)")+
  ggtitle("B")
p4 <- qualplot+
  #  ylab("PCO2 (12.3%)")+
  geom_point(aes(col=qualplate))
p6 <- qualplot+
  geom_point(aes(col=pop(qualgenind)))
p5 <- qualplot+
  #  xlab("PCO1 (10.7%)")+
  geom_point(aes(col=qualmiss))+
  geom_text(data=qualpco$li[qualsub,], aes(label=c("Lssn_M11","BMT_M05","NDunes_F03")))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p4,p5,p6, nrow=1)




#### 5A: Minimum quality-filtered w/ outliers removed ####
# extract tab
nooutgenind <- vcfR2genind(vcf5a)
colnames(vcf5a@gt)[!colnames(vcf5a@gt)%in%indNames(nooutgenind)]#MACAtail, MACA_HU_F05
noouttab <- tab(nooutgenind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(nooutgenind) <- popmap$pop[match(indNames(nooutgenind), popmap$ID)]
nooutplate <- c(rep("Plate1",37),rep("Plate2",26),rep("Plate3",29))
nooutmiss <- apply(noouttab, 1, function(x) sum(is.na(x))/length(x))

# run pco
nooutpco <- dudi.pco(dist(noouttab), scannf=FALSE, nf=50)
summary(nooutpco)

# plots
nooutplot <- ggplot(nooutpco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (14.6%)")+
  ylab("PCO2 (10.2%)")+
  ggtitle("C")
p7 <- nooutplot+
  #  ylab("PCO2 (10.2%)")+
  geom_point(aes(col=nooutplate))
p9 <- nooutplot+
  geom_point(aes(col=pop(nooutgenind)))
p8 <- nooutplot+
  #  xlab("PCO1 (14.6%)")+
  geom_point(aes(col=nooutmiss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p7,p8,p9, nrow=1)



#### 5B: Quality, no outliers, 30% genotyping ####
# extract tab
geno30genind <- vcfR2genind(vcf5b)
colnames(vcf5b@gt)[!colnames(vcf5b@gt)%in%indNames(geno30genind)]#MA-16-BMT8, MACAtail, MACA_HU_F05
geno30tab <- tab(geno30genind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(geno30genind) <- popmap$pop[match(indNames(geno30genind), paste0("0_",popmap$ID))]
geno30plate <- c(rep("Plate1",36),rep("Plate2",26),rep("Plate3",29))
geno30miss <- apply(geno30tab, 1, function(x) sum(is.na(x))/length(x))

# run pco
geno30pco <- dudi.pco(dist(geno30tab), scannf=FALSE, nf=50)
summary(geno30pco)

# plots
geno30plot <- ggplot(geno30pco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (15.8%)")+
  ylab("PCO2 (7.5%)")+
  ggtitle("D")
p10 <- geno30plot+
  #  ylab("PCO2 (7.5%)")+
  geom_point(aes(col=geno30plate))
p12 <- geno30plot+
  geom_point(aes(col=pop(geno30genind)))
p11 <- geno30plot+
  #  xlab("PCO1 (15.8%)")+
  geom_point(aes(col=geno30miss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p10,p11,p12, nrow=1)



#### 5C: Quality, no outliers, 30% genotyping, 50% missingness per ind ####
# extract tab
mind50genind <- vcfR2genind(vcf5c)
colnames(vcf5c@gt)[!colnames(vcf5c@gt)%in%indNames(mind50genind)]#
mind50tab <- tab(mind50genind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(mind50genind) <- popmap$pop[match(indNames(mind50genind), popmap$ID)]
pop(mind50genind)[c(36:38)] <- "SDunes"
mind50plate <- c(rep("Plate1",28),rep("Plate2",23),rep("Plate3",24))
mind50miss <- apply(mind50tab, 1, function(x) sum(is.na(x))/length(x))

# run pco
mind50pco <- dudi.pco(dist(mind50tab), scannf=FALSE, nf=50)
summary(mind50pco)

# plots
mind50plot <- ggplot(mind50pco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (16.5%)")+
  ylab("PCO2 (7.8%)")+
  ggtitle("E")
p13 <- mind50plot+
  #  ylab("PCO2 (7.8%)")+
  geom_point(aes(col=mind50plate))
p15 <- mind50plot+
  geom_point(aes(col=pop(mind50genind)))+
  scale_color_manual(values=hue_pal()(11)[match(levels(pop(mind50genind)), levels(pop(rawgenind)))])
p14 <- mind50plot+
  #  xlab("PCO1 (16.5%)")+
  geom_point(aes(col=mind50miss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p13,p14,p15, nrow=1)



#### 10: Fully Filtered ####
# extract tab
filtgenind <- vcfR2genind(vcf10)
colnames(vcf10@gt)[!colnames(vcf10@gt)%in%indNames(filtgenind)]
filttab <- tab(filtgenind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(filtgenind) <- popmap$pop[match(indNames(filtgenind), paste(popmap$ID,popmap$ID,sep="_"))]
filtplate <- c(rep("Plate1",26),rep("Plate2",16),rep("Plate3",21))
filtmiss <- apply(filttab, 1, function(x) sum(is.na(x))/length(x))

# run pco
filtpco <- dudi.pco(dist(filttab), scannf=FALSE, nf=50)
summary(filtpco)

# plots
filtplot <- ggplot(filtpco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (16.9%)")+
  ylab("PCO2 (7.9%)")+
  ggtitle("F")
p16 <- filtplot+
  #  ylab("PCO2 (7.9%)")+
  geom_point(aes(col=filtplate))
p18 <- filtplot+
  geom_point(aes(col=pop(filtgenind)))+
  scale_color_manual(values=hue_pal()(11)[match(levels(pop(filtgenind)), levels(pop(rawgenind)))])
p17 <- filtplot+
  #  xlab("PCO1 (16.9%)")+
  geom_point(aes(col=filtmiss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p16,p17,p18, nrow=1)



#### 11: Low missingness dataset (5% geno, 10% mind) ####
# extract tab
lowmissgenind <- vcfR2genind(vcf12)
colnames(vcf12@gt)[!colnames(vcf12@gt)%in%indNames(lowmissgenind)]#MACA_HU_F05A, MACAtail
lowmisstab <- tab(lowmissgenind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(lowmissgenind) <- popmap$pop[match(indNames(lowmissgenind), paste(popmap$ID,popmap$ID,sep="_"))]
pop(lowmissgenind)[c(4,27,44,58)] <- c("BlueMountains","SDunes","SouthernOregon","SDunes")
lowmissplate <- c(rep("Plate1",22),rep("Plate2",16),rep("Plate3",20))
lowmissmiss <- apply(lowmisstab, 1, function(x) sum(is.na(x))/length(x))

# run pco
lowmisspco <- dudi.pco(dist(lowmisstab), scannf=FALSE, nf=50)
summary(lowmisspco)

# plots
lowmissplot <- ggplot(lowmisspco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (19.3%)")+
  ylab("PCO2 (7.4%)")+
  ggtitle("G")
p19 <- lowmissplot+
  #  ylab("PCO2 (7.4%)")+
  geom_point(aes(col=lowmissplate))
p21 <- lowmissplot+
  geom_point(aes(col=pop(lowmissgenind)))+
  scale_color_manual(values=hue_pal()(11)[match(levels(pop(lowmissgenind)), levels(pop(rawgenind)))])
p20 <- lowmissplot+
  #  xlab("PCO1 (19.3%)")+
  geom_point(aes(col=lowmissmiss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p19,p20,p21, nrow=1)



#### 12: High missingness +  LD-pruning ####
# extract tab
ld10genind <- vcfR2genind(vcf11)
colnames(vcf11@gt)[!colnames(vcf11@gt)%in%indNames(ld10genind)]#MACA_HU_F05A, MACAtail
ld10tab <- tab(ld10genind, freq=T, NA.method="asis")

# grab ids, plates, missingness
pop(ld10genind) <- popmap$pop[match(indNames(ld10genind), paste("0",popmap$ID,popmap$ID,sep="_"))]
ld10plate <- c(rep("Plate1",26),rep("Plate2",16),rep("Plate3",21))
ld10miss <- apply(ld10tab, 1, function(x) sum(is.na(x))/length(x))

# run pco
ld10pco <- dudi.pco(dist(ld10tab), scannf=FALSE, nf=50)
summary(ld10pco)

# plots
ld10plot <- ggplot(ld10pco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (17.1%)")+
  ylab("PCO2 (6.5%)")+
  ggtitle("H")
p22 <- ld10plot+
  #  ylab("PCO2 (6.5%)")+
  geom_point(aes(col=ld10plate))
p24 <- ld10plot+
  geom_point(aes(col=pop(ld10genind)))+
  scale_color_manual(values=hue_pal()(11)[match(levels(pop(lowmissgenind)), levels(pop(rawgenind)))])
p23 <- ld10plot+
  #  xlab("PCO1 (17.1%)")+
  geom_point(aes(col=ld10miss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p22,p23,p24, nrow=1)



#### 12B: Low missingness + LD-pruning ####
# extract tab
ld11genind <- vcfR2genind(vcf12b)
colnames(vcf12b@gt)[!colnames(vcf12b@gt)%in%indNames(ld11genind)]#MACA_HU_F05A, MACAtail
ld11tab <- tab(ld11genind, freq=T, NA.method="mean")

# grab ids, plates, missingness
pop(ld11genind) <- pop(lowmissgenind)
ld11plate <- lowmissplate
ld11miss <- apply(ld11tab, 1, function(x) sum(is.na(x))/length(x))

# run pco
ld11pco <- dudi.pco(dist(ld11tab), scannf=FALSE, nf=50)
summary(ld11pco)

# plots
ld11plot <- ggplot(ld11pco$li, aes(x=A1, y=A2, size=2, alpha=0.7))+
  theme_bw()+
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab("PCO1 (19.0%)")+
  ylab("PCO2 (7.0%)")+
  ggtitle("I")
p25 <- ld11plot+
  #  ylab("PCO2 (7.0%)")+
  geom_point(aes(col=ld11plate))
p27 <- ld11plot+
  geom_point(aes(col=pop(ld11genind)))+
  scale_color_manual(values=hue_pal()(11)[match(levels(pop(ld11genind)), levels(pop(rawgenind)))])
p26 <- ld11plot+
  #  xlab("PCO1 (19.0%)")+
  geom_point(aes(col=ld11miss))+
  theme(legend.position="right")+
  scale_color_continuous(limits=c(0,1))+
  guides(size="none", alpha="none")
ggpubr::ggarrange(p25,p26,p27, nrow=1)



###############################
#### FINAL FILTERING PLOTS ####
###############################
leg <- ggpubr::get_legend(p1+
                            theme(legend.position="right")+
                            guides(size="none", alpha="none"))
ggpubr::ggarrange(p1,p4,p7,p10,p13,p16,p19,p22,p25, legend.grob=leg, legend="right")

ggpubr::ggarrange(p2,p5,p8,p11,p14,p17,p20,p23,p26, common.legend=TRUE, legend="right")

leg <- ggpubr::get_legend(p3+
                            theme(legend.position="right")+
                            guides(size="none",alpha="none"))
ggpubr::ggarrange(p3,p6,p9,p12,p15,p18,p21,p24,p27, legend.grob=leg, legend="right")

