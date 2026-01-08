# load dependencies
library(vcfR)
library(hierfstat)
library(SNPRelate)
library(corrplot)
library(PopGenReport)
library(GenoPop)
library(stringr)


## Read in data and convert to genind
vcf11 <- vcfR::read.vcfR("LDprunedRAD.vcf")
genind11 <- vcfR::vcfR2genind(vcf11)

# assign populations based on individual IDs
adegenet::pop(genind11) <- stringr::str_split(adegenet::indNames(genind11), "_", n=2, simplify=TRUE)[,1]
table(pop(genind11))#check
# factor pops for ease...
pop(genind11)[pop(genind11)=="Border"] <- "NCali" #reassign one border individual to NCALI
pop(genind11) <- factor(pop(genind11), levels=c("NDunes","SDunes","SOreg","NCali","MapleCrk",
                                                "ORcasc","Lssn","WAcasc","ORbmt","WY","CO"))

# designate lineage for each individuals
lineage <- as.character(pop(genind11))
lineage[lineage%in%c("NDunes","SDunes","SOreg","NCali","MapleCrk","Border")] <-"Humboldt"
lineage[lineage%in%c("ORcasc","Lssn")] <- "WMontane"
lineage[lineage%in%c("WAcasc","CO","ORbmt","WY")] <- "EMontane"


## Subset for Humboldt martens - include monomorphic SNPs!
# subset VCF to Humboldt martens
vcfcoma11 <- vcf11[,c(1,which(lineage=="Humboldt")+1)]
# calc MAFs
maf <- vcfR::maf(vcfcoma11)
# subset to MAFs>0
vcfcoma11 <- vcfcoma11[maf[,4]>0,]; vcfcoma11
# convert to genind
genindcoma11 <- vcfR::vcfR2genind(vcfcoma11)
# designate pops
adegenet::pop(genindcoma11) <- stringr::str_split(adegenet::indNames(genindcoma11), "_", n=2, simplify=TRUE)[,1]
genindcoma11#39 x 9239 SNPs remaining


## create subset for eastern and western montane martens following same process
vcfmont11 <- vcf11[,c(1,which(lineage=="EMontane")+1)]
maf <- vcfR::maf(vcfmont11)
vcfmont11 <- vcfmont11[maf[,4]>0,]
genindemont <- vcfR2genind(vcfmont11)
adegenet::pop(genindemont) <- stringr::str_split(adegenet::indNames(genindemont), "_", n=2, simplify=TRUE)[,1]
genindemont#10168 SNPs, 16 individuals

vcfmont11 <- vcf11[,c(1,which(lineage=="WMontane")+1)]
maf <- vcfR::maf(vcfmont11)
vcfmont11 <- vcfmont11[maf[,4]>0,]
genindwmont <- vcfR2genind(vcfmont11)
adegenet::pop(genindwmont) <- stringr::str_split(adegenet::indNames(genindwmont), "_", n=2, simplify=TRUE)[,1]
genindwmont#7821 SNPs, 8 individuals



#### POPULATION-LEVEL FST, FIS ####
# Fst for each lineage relative to total
macafs <- fs.dosage(genind11, lineage)
plot(macafs)
macafs$Fs

# Fsts relative to all Humboldt...
fs.dosage(genindcoma11, pop(genindcoma11))

# Fsts relative to all montane....
fs.dosage(genindemont, pop(genindemont))
fs.dosage(genindwmont, pop(genindwmont))



#### PAIRWISE FST ####
# calculate Fsts and 95% bootstrap C.I.s
fst <- hierfstat::pairwise.neifst(genind11)
fst_999bs <- hierfstat::boot.ppbetas(genind11, nboot=999)

# plot
corrplot.mixed(fst, lower="number", upper="shade", upper.col=rev(COL2('RdBu',20)), 
               lower.col='black', col.lim=c(0,1), diag='n')



#### OBSERVED AND EXPECTED HETEROZYGOSITY (Hs, Ho) ####
# all pops
hierfstat::Ho(genind11)
hierfstat::Hs(genind11)

# lineages
genind_lins <- genind11
pop(genind_lins) <- lineage
hierfstat::Ho(genind_lins)
hierfstat::Hs(genind_lins)



#### ALLELIC RICHNESS, PRIVATE ALLELES, & PROP SNPs SHARED (Ar) ####
## calculate overall allelic richness per pop
ar <- hierfstat::allelic.richness(genind11, min.n=3)
summary(ar$Ar)
boxplot(ar$Ar)

# and by lineage
ar2 <- hierfstat::allelic.richness(genind_lins)
summary(ar2$Ar)



#### PRIVATE ALLELIC RICHNESS PER POP (Ap) ####
# extract private alleles from genind
ad <- PopGenReport::allele.dist(genind11)
pvt <- as.data.frame(t(as.data.frame(ad$private.alleles)))
pvt <- pvt[which(!is.na(pvt$Population)),]
pvt$Allele <- as.numeric(pvt$Allele)
aggregate(Allele~Population, data=pvt, FUN=function(x) sum(x,na.rm=T))
genindIDs <- unlist(lapply(colnames(genind11@tab), function(x) strsplit(x,":")[[1]][1]))
pvtIDs <- unlist(lapply(rownames(pvt), function(x) gsub("X","", strsplit(x,"\\.")[[1]][1])))
pvtsub <- genind11[,genindIDs%in%pvtIDs]; pvtsub; nrow(pvt)

# run allelic richness for private alleles only
pa <- hierfstat::allelic.richness(pvtsub, min.n=3)
summary(pa$Ar)
pop(pvtsub) <- lineage
pa2 <- hierfstat::allelic.richness(pvtsub)
summary(pa2$Ar)

## shared alleles per population
nb <- hierfstat::nb.alleles(genind11); nb



#### AUTOSOMAL HETEROZYGOSITY (Ha) ####
tab <- tab(genind11, freq=T, NA.method="asis")
loci <- read.csv("populations.loci.csv") #this is the output of populations --fasta-loci converted to a 2-column CSV (Name, Sequence)
loci$Length <- unlist(lapply(loci$Sequence, function(x) nchar(x)))

# subset loci to those included in VCF
fixes <- as.data.frame(vcf11@fix)
fixes$CHROM <- unlist(lapply(fixes$ID, function(x) strsplit(x,":")[[1]][1]))
chroms <- paste0("CLocus_", fixes$CHROM)
locisub <- loci[loci$Locus%in%chroms,]

# extract # HET per individual
numhet <- apply(tab, 1, function(x) sum(x==0.5,na.rm=T))

# calc Ha per individual
ha_df <- data.frame(ID=rownames(tab), Ha=NA)
for (ind in rownames(tab)){
  # subset to genotyped loci
  l <- colnames(tab)[!is.na(tab[ind,])]
  lchrom <- unlist(lapply(l, function(x) gsub("X","",strsplit(x,":")[[1]][1])))
  lsub <- locisub[locisub$Locus %in% paste0("CLocus_",lchrom),]
  # calculate Ha: Het SNPs / total BP- this is an overestimate because only includes variant sites!!!
  h <- numhet[ind]
  a <- sum(lsub$Length)
  ha_df$Ha[ha_df$ID==ind] <- h / a
}#for

# means per pop
ha_df$Pop <- pop(genind11)
aggregate(Ha~Pop, ha_df, FUN=mean)
aggregate(ha_df$Ha~lineage, FUN=mean)
## NOTE: THESE ARE OVERESTIMATES B/C THEY DON'T INCLUDE NONVARIANT SITES!!!



#### INDIVIDUAL Hsnp ####
## SNP-Based Heterozygosity
numgeno <- apply(tab, 1, function(x) sum(!is.na(x)))
hsnp <- numhet / numgeno
hsnp <- data.frame(ID=indNames(genind11),
                   pop=pop(genind11),
                   HSNP=hsnp)

# mean per population
aggregate(HSNP~pop, hsnp, mean)
aggregate(hsnp$HSNP~lineage, FUN=mean)



#### INDIVIDUAL INBREEDING COEFFICIENTS ####
# Weir and Goudet beta estimator - this is equivalent to plink FHOM estimator
betas <- hierfstat::beta.dosage(genind11, inb=TRUE)#Mb=TRUE  to report mean matching
hsnp$Fbeta <- as.vector(diag(betas))

# Fis based on Weir and Goudet
aggregate(Fbeta~pop, hsnp, FUN=mean)
aggregate(hsnp$Fbeta~lineage, FUN=mean)

# check correlation 
cor(hsnp$HSNP, hsnp$Fbeta)#-1



#### HSNP AND FBETA PLOTS ####
# rename pops for plotting
het <- hsnp
het$pop <- as.character(het$pop)
het$pop[het$pop=="NDunes"] <- "N Dunes"
het$pop[het$pop=="SDunes"] <- "S Dunes"
het$pop[het$pop=="SOreg"] <- "S Oregon"
het$pop[het$pop=="NCali"] <- "N California"
het$pop[het$pop=="MapleCrk"] <- "Maple Crk"
het$pop[het$pop=="CO"] <- "N Colorado"
het$pop[het$pop=="WY"] <- "Wyoming"
het$pop[het$pop=="ORbmt"] <- "Blue Mtns"
het$pop[het$pop=="WAcasc"] <- "WA Cascades"
het$pop[het$pop=="Lssn"] <- "Lassen"
het$pop[het$pop=="ORcasc"] <- "OR Cascades"
het$pop <- factor(het$pop, levels=c("N Colorado","Wyoming","Blue Mtns","WA Cascades","Lassen","OR Cascades",
                  "Maple Crk","N California","S Oregon", "S Dunes", "N Dunes","Humboldt"))

# add lineage-level to DF so that these can also be plotted
hetsub <- het
hetsub$pop[hetsub$pop%in%c("Maple Crk","N California","S Oregon", "S Dunes", "N Dunes")] <- "Humboldt"
hetsub$pop[hetsub$pop%in%c("Lassen","OR Cascades")] <- "W Montane"
hetsub$pop[hetsub$pop%in%c("Blue Mtns","WA Cascades","Wyoming","N Colorado")] <- "E Montane"
hetall <- rbind(hetsub, het)

# factor populations/lineages for plotting order
hetall$pop <- factor(hetall$pop, levels=c("N Colorado","Wyoming","Blue Mtns","WA Cascades","E Montane","Lassen","OR Cascades","W Montane",
                                          "Maple Crk","N California","S Oregon", "S Dunes", "N Dunes","Humboldt"))

# set up plotting colors
# color palette
cols <- c(#"#000000",#black-UNK
          "#db6d00",#orange-CO
          "gray20",#"#252525",#darkgray-WY
          "gray60",#676767",#lightgray-ORbmt
          "gray90",#"#ffffff",#white-WAcasc
          "red",#east montane
          "mediumpurple4",#171723",#navy-Lssn
          "darkcyan",#"#004949",#dark teal-ORcasc
          "mediumpurple4",#west montane
          #"aquamarine1",#009999",#light teal-Applegate
          "#22cf22",#bright green-NCali
          "darkviolet",#"#490092",#purple-MapleCrk
          #"steelblue1",#"#006ddb",#sky blue-Border
          #"#b66dff",#light purple-SCoast
          "#ff6db6",#bright pink-SOreg
          "#920000",#dark red-SDunes
          "#ffdf4d",#yellow-NDunes
          #,"#8f4e00")#brown-NCoast
          "blue")#humboldt

# # boxplots
# ggplot(hetall)+
#   geom_boxplot(aes(x=pop, y=HSNP, fill=pop))+
#   theme_classic()+
#   scale_fill_manual(values=cols)+
#   theme(legend.position="none", text=element_text(size=16))+
#   xlab(element_blank())+
#   ylab(expression(H[SNP]))+
#   coord_flip()
# ggplot(hetall)+
#   geom_boxplot(aes(x=pop, y=Fbeta, fill=pop))+
#   theme_classic()+
#   scale_fill_manual(values=cols)+
#   theme(legend.position="none", text=element_text(size=16))+
#   xlab(element_blank())+
#   ylab(expression(hat(F[beta])))+
#   coord_flip()

## plot raw genetic diversity estimates per population
ggplot(het)+
  geom_point(aes(x=pop, y=HSNP, col=pop), size=5, alpha=0.7)+
  scale_color_manual(values=cols[-c(5,8,14)])+
  theme_bw()+
  coord_flip()+
  xlab(element_blank())+
  ylab(expression(H[SNP]))+
  theme(legend.position="none", text=element_text(size=16))

## plot individual genetic diversity x inbreeding coefficients by lineage
ggplot(hetall[which(hetall$pop%in%c("Humboldt","E Montane","W Montane")),])+
  geom_point(aes(x=HSNP, y=Fbeta, col=pop), size=4, alpha=0.7)+
  scale_color_manual(values=c("dodgerblue","purple4","red"),
                     breaks=c("Humboldt","E Montane","W Montane"))+
  theme_bw()+
  #coord_flip()+
  ylab(expression(F[beta]))+
  xlab(expression(H[SNP]))+
  theme(legend.position="inside", legend.position.inside=c(0.8,0.8), text=element_text(size=16), legend.title=element_blank())

# boxplots for Hsnp by lineage w/ sample sizes in xlab
N <- as.data.frame(table(hetall$pop))
hetall$N <- paste("N=",N$Freq[match(hetall$pop, N$Var1)],sep="")
ggplot(hetall[which(hetall$pop%in%c("Humboldt","E Montane","W Montane")),])+
  geom_boxplot(aes(x=N, y=HSNP, fill=pop))+
  scale_fill_manual(values=c("blue","purple","red"), breaks=c("Humboldt","W Montane","E Montane"))+
  theme_classic()+
  xlab(element_blank())+
  ylab(expression(H[SNP]))+
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=45,vjust=0.5),
        text=element_text(size=12), legend.position="inside", legend.position.inside=c(0.8,0.8),
        legend.background=element_rect(fill=NULL, colour="black"),legend.title=element_blank())



##### CHECK HOW MANY LOCI HAVE >2 SNPS.... ####
tbl <- table(chroms)
tbl <- as.data.frame(tbl)
hist(tbl$Freq)
sum(tbl$Freq>2)/12389 #<5%



#### HETEROZYGOSITY VS DEPTH ####
## included in SI
# calculate mean depths per sample
dp <- vcfR::extract.gt(vcf11, element="DP")
sample_depths <- apply(dp, 2, FUN=function(x) mean(as.numeric(x),na.rm=T))
hist(sample_depths,main="Sequencing Depth per Sample",xlab="Mean SNP depth")
hetall$MeanDP <- sample_depths[match(names(sample_depths), hetall$ID)]

# add plate/batch to DF
plate <- c(rep("Plate1",26), rep("Plate2",16), rep("Plate3",21))
hetall$Plate <- plate[match(indNames(genind11), hetall$ID)]

# plot mean depth vs. heterozygosity to check for effects of low depth on gendiv estimate
which(hetall$HSNP>0.25)
ggplot(hetall[1:63,], aes(x=MeanDP, y=HSNP))+
   geom_point(size=2, aes(col=Plate))+
   geom_text(data=hetall[c(5,11,16),], aes(y=HSNP+0.01, label=ID), position=position_jitter())+
   theme_classic()+
   ylab(expression(H[SNP]))+
   xlab("Sequencing depth")



#### MANTEL'S TEST ####
library(sp)
library(geosphere)
library(vegan)
# read in sample coordinates
popmap <- read.csv("../SampleLocations/Marten_Hets_Fhom_coords.csv", stringsAsFactors=FALSE, fileEncoding='Latin1', header=T)
# find matching rows between RADseq data and spatial CSV
matches <- match(indNames(genind11), paste(popmap$MergeIDraw, popmap$MergeIDraw, sep="_"))
indNames(genind11)[is.na(matches)]#check which are missing

# extract coordinates
gtcoords <- data.frame(Longitude=popmap$Longitude[match(indNames(finalgenind), popmap$MergeID)],
                       Latitude=popmap$Latitude[match(indNames(finalgenind), popmap$MergeID)])
gtcoords$Longitude <- as.numeric(gtcoords$Longitude)
sum(is.na(gtcoords$Latitude))# check for missing

## convert lat/long coords to a projection that preserves distance
crs <- CRS("EPSG:4326")
sp <- sp::SpatialPointsDataFrame(gtcoords[,c("Longitude","Latitude")], gtcoords, proj4string=crs)
sp <- spTransform(sp, CRS("ESRI:102005"))#USA contiguous equidistant conic
gtcoords <- as.data.frame(coordinates(sp))
names(gtcoords) <- c("x","y")

## Run Mantel's test
# grab centroid coordinate for each population/region
gtcoords$pop <- pop(genind11)
popcoords <- aggregate(data=gtcoords, Longitude~pop, FUN=mean)
popcoords$Latitude <- aggregate(data=gtcoords, Latitude~pop, FUN=mean)[,2]
# calculate pairwise distances (meters) between all populations
popdist <- geosphere::distm(popcoords[,c("Longitude","Latitude")])
# run mantels test 
vegan::mantel(fst, popdist, method="pearson")
#mc <- vegan::mantel.correlog(fst, popdist, cutoff=TRUE)
#plot(mc)
#mc
# ggplot()+
#   geom_point(aes(x=as.vector(popdist)/1000, 
#                  y=as.vector(fst)))+ 
#   #col=rep(popcoords$pop, nrow(popcoords))))+
#   xlab("Geographic distance")+
#   ylab("Fst")+
#   theme_classic()+
#   geom_vline(aes(xintercept=140), col="red")



#### PLOT PCO AXES OVER SPACE ####
library(akima)
library(maps)
library(sf)
# run pco
pco <- dudi.pco(dist(tab(genind11)), scannf=FALSE, nf=2)

# interpolate PCO1 over a grid
## note: a couple samples didn't have precise coordinates and are excluded here
x <- gtcoords$Longitude
y <- gtcoords$Latitude

makeGrid <- function(x){
  # grab min/max values
  min <- min(x)
  max <- max(x)
  range <- max-min
  # add 10% to either end
  min <- min-(range/10)
  max <- max+(range/10)
  # divide range into 50 evenly spaced points
  grid <- seq(min, max, length.out=200)
  return(grid)
}#makeGrid

# load in basemap of state boundaries & transform proj
states <- maps::map_data("state",)
statesp <- sp::SpatialPointsDataFrame(states[,c("long","lat")], states, proj4string=CRS("EPSG:4326"))
statesp <- sp::spTransform(statesp, CRS("ESRI:102005"))
states$long <- coordinates(statesp)[,1]
states$lat <- coordinates(statesp)[,2]

# interpolate PCO1 over grid
nas <- which(is.na(x))
x <- x[-nas]
y <- y[-nas]
pco1 <- pco$li$A1[-nas]
img <- interp(x, y, pco1, duplicate="median", xo=makeGrid(x), yo=makeGrid(y))

ggplot()+
  geom_map(data=states, map=states, aes(map_id=region), fill=NA, col='black')+
  geom_contour_filled(aes(x=rep(img$x, nrow(img$z)), 
                          y=rep(img$y, each=ncol(img$z)),
                          z=as.numeric(img$z)), alpha=0.5, show.legend=FALSE, bins=10)+
  geom_point(aes(x=x, y=y, col=pco1), alpha=0.8, size=3)+
  theme_bw()+
  #xlab("Longitude")+
  #ylab("Latitude")+
  labs(col="PCO1")+
  #ggtitle("PCoA Axis1 in Geographic Space")+
  coord_fixed()+
  scale_color_viridis()+
  scale_fill_viridis(discrete=T)+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())+
  coord_sf(crs=CRS("ESRI:102005"),xlim=c(-2350546.916,-840711.114), ylim=c(200000,1352610.596))


## rerun with Humboldt martens only...
coma_pco <- dudi.pco(dist(tab(genindcoma11)), scannf=FALSE, nf=2)
gtcoords_coma <- gtcoords[pop(genind11)%in%c("SDunes","MapleCrk","NCali","NDunes","SOreg"),]
x <- gtcoords_coma$Longitude; y <- gtcoords_coma$Latitude
coma_img <- akima::interp(x[-which(is.na(x))], y[-which(is.na(x))], coma_pco$li$A1[-which(is.na(x))], duplicate="median")




#### Checking how many SNPs are variable across Humboldt marten pops ####
vcf <- read.vcfR("../RADseq_Plate3_Results/VCFs/11_LD50560.vcf")

grabSNPs <- function(vcf, pop, threshold, gi=genind){
  sub_vcf <- vcf[,c(1, which(pop(gi)==pop)+1)]
  mafs <- as.data.frame(maf(sub_vcf))
  sub_vcf <- sub_vcf[which(mafs$Frequency>threshold & mafs$Frequency<(1-threshold)),]
  fixes <- as.data.frame(sub_vcf@fix)
  return(fixes$ID)
}#grabSNPs

ndunes_snps <- grabSNPs(vcf, "NDunes", 0.2); length(ndunes_snps)#4390 SNPs, 2638>0.2
sdunes_snps <- grabSNPs(vcf, "SDunes", 0.2); length(sdunes_snps)#5949, 2871>0.2
soreg_snps <- grabSNPs(vcf, "SOreg", 0.2); length(soreg_snps)#6619, 3271>0.2
ncali_snps <- grabSNPs(vcf, "NCali", 0.2); length(ncali_snps)#6736, 3144>0.2

tbl <- table(c(ndunes_snps, sdunes_snps, soreg_snps, ncali_snps))
sum(tbl==4)#2605 total are polymorphic across ALL sites
