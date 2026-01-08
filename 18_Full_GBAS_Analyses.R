# load dependencies
library(vcfR)
library(adegenet)
library(ggplot2)
library(ggpubr)
library(sp)
library(maps)
library(akima)
library(viridis)
library(RColorBrewer)
library(popgraph)
library(phangorn)
library(ape)
library(poppr)
library(magrittr)
library(ggtree)
library(colorspace)


################################
##### READ IN GENOTYPE DATA ####
################################
# load genotypes
inds <- read.csv("data/GBAS_sample_identification.csv")
#subset to individuals selected to use & with plotID
inds <- inds[which(inds$Use==1),]; nrow(inds)#150
#remove sex loci
inds <- inds[,-which(names(inds)%in%c("SRY_MAAM.1","SRY_MAAM.4"))] 

# convert to plink format
{
  # subset to SNPs and clean up NAs
  genos <- inds[,c(13:113)]
  genos[genos=="0"] <- NA
  genos[genos=="NA"] <- NA
  genos[genos=="-"] <- NA
  genos[genos=="FALSE"] <- NA
  sort(unique(unlist(genos)))#check
  
  # convert to plink format
  genos[genos=="AA"] <- "A A"
  genos[genos=="AC"] <- "A C"
  genos[genos=="AG"] <- "A G"
  genos[genos=="AT"] <- "A T"
  genos[genos=="CA"] <- "C A"
  genos[genos=="CC"] <- "C C"
  genos[genos=="CG"] <- "C G"
  genos[genos=="CT"] <- "C T"
  genos[genos=="GA"] <- "G A"
  genos[genos=="GC"] <- "G C"
  genos[genos=="GG"] <- "G G"
  genos[genos=="GT"] <- "G T"
  genos[genos=="TA"] <- "T A"
  genos[genos=="TC"] <- "T C"
  genos[genos=="TG"] <- "T G"
  genos[genos=="TT"] <- "T T"
  genos[is.na(genos)] <- "0 0"
  
  # copy into DF
  dim(genos)# check sample size vs. # SNPs
  inds[,13:113] <- genos
  
  # clean up locus names
  names(inds) <- gsub("\\.", "_",names(inds))
  inds$Pop[is.na(inds$Pop)]#check
  #View(plink)
  #write.csv(plink, "individualIDs_17jun2025_plink.csv")
}

## convert PLINK format to genind
# remove microhaplotypes that are fully linked
inds <- inds[,-which(names(inds)%in%c("MACA_12025_3_158","MACA_131823_2_52","MACA_84320_4_72"))]

# remove samples that say not to use (this is using full dataset since some RAD samples failed to genotype at 75SNP threshold)
inds <- inds[which(inds$Use==1),]; nrow(inds)#150
inds$SampleID[inds$PlotID%in%c(NA,"")]

# convert to genind
names(inds)
snpgenind <- df2genind(inds[,13:110], sep=" ", NA.char="0 0", ind.names=inds$PlotID, pop=inds$Pop)
snpgenind#150 martens x 98 loci
unique(pop(snpgenind))#check
pop(snpgenind)#check
indNames(snpgenind)#check

# set up plotting variables
{
  pops <- as.character(pop(snpgenind))
  pops[pops=="Blue Mountains"] <- "Blue Mtns"
  pops[pops=="S Oregon?"] <- "S Oregon"
  pops[pops=="North Coast"] <- "N Coast OR"
  pops[pops%in%c("S Coast")] <- "S Coast OR"
  pops[pops=="Border"] <- "OR/CA Border"
  pops <- factor(pops, levels=c("N Coast OR","N Dunes", "S Dunes", 
                                "S Coast OR", "S Oregon", "OR/CA Border", 
                                "Trinidad", "N California", "Ashland", "OR Cascades", "WA Cascades", 
                                "Lassen", "Blue Mtns",
                                "Wyoming", "Colorado", "Unknown","MAAM")) 
  pops
}

# pull sample ID from plotID
type <- unlist(lapply(indNames(snpgenind), function(x) strsplit(x,"\\.")[[1]][2]))
type[type%in%c("rk","tr")] <- "tissue"
type[type%in%c("sc","hr")] <- "scat"
indNames(snpgenind)[is.na(type)]
type[is.na(type)] <- "tissue"



##################
#### RUN PCoA ####
##################
{
  dist <- dist(tab(snpgenind[-which(pops=="MAAM"),], freq=T, method="as.is"))#remove americana
  pco <- dudi.pco(dist, nf=200, scannf=FALSE)
  #s.class(-pco$li, fac=pop(snpgenind)[-which(pops=="MAAM")], col=rainbow(15), cstar=0, axesell=0, label=NULL)
  #add.scatter.eig(pco$eig, posi="bottomright", xax=1, yax=2)
  #plot(-40:0, -40:0, pch=0)
  #legend(x=-40, y=1, legend=levels(pop(genind_sub)), col=rainbow(15), pch=16)
  summary(pco)
  dim(pco$li)

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
  pops <- factor(pops, levels=c("N Coast OR","N Dunes","S Dunes","S Oregon","S Coast OR","OR/CA Border","Trinidad",
                                "N California","Ashland","OR Cascades","Lassen","WA Cascades","Blue Mtns","Wyoming","Colorado","Unknown","MAAM"))
  pops[pops=="Unknown"] <- "S Dunes"
  ggplot(data=pco$li, aes(x=A1, y=A2, col=pops[-which(pops=="MAAM")]))+
    geom_hline(aes(yintercept=0))+
    geom_vline(aes(xintercept=0))+
    geom_point(aes(pch=type[-which(pops=="MAAM")]), size=5, alpha=0.8)+
    #stat_chull(aes(fill=pops[-which(pops=="MAAM")]), geom="polygon",lwd=0,alpha=0.2)+
    theme_bw()+
    #ggtitle("SNP Panel (98 SNPs, 148 martens)")+
    scale_color_manual(breaks=levels(pops[-which(pops=="MAAM")]), values=rev(pal))+#values=c(rainbow(15),"grey"))+
    theme_classic()+
    xlab("PCO1 (15.1%)")+
    ylab("PCO2 (12.6%)")+
    labs(col="Population")+
    ggtitle("GBAS (N=148, SNPs=98)")+
    guides(fill="none", shape="none")+
    theme(text=element_text(size=20), axis.text=element_blank(), axis.ticks=element_blank(),
          plot.title=element_text(face="bold",size=20), plot.background=element_rect(fill=NA, color='black'), 
          axis.line=element_blank(), legend.position="none")#panel.grid.major=element_line(color='gray70',linewidth=0.5), 
    #coord_fixed()
  #add.scatter.eig(pco$eig, posi="topright", xax=1, yax=2, ratio=0.3)
}
#ggplot(snpcomapco$li)+geom_point(aes(x=A1,y=A2,col=pop(snpcoma)), size=2,alpha=0.8)+theme_bw()+coord_fixed()+xlab("PCO1 (19.3%)")+ylab("PCO2 (9.7%)")+labs(color="Population")+scale_color_manual(values=rev(cols), breaks=levels(df$Pop1))+theme(axis.text=element_blank(),axis.ticks = element_blank(), text=element_text(size=18))



#### PCoA - HUMBOLDT MARTENS ONLY #### 
{
  humb <- c("N Coast OR","N Dunes", "S Dunes", "S Coast OR", "S Oregon", "OR/CA Border", 
            "Trinidad", "N California", "Ashland","Unknown")
  dist <- dist(tab(snpgenind[which(pops%in%humb),], freq=T, method="as.is"))#remove americana
  pco <- dudi.pco(dist, nf=200, scannf=FALSE)
  summary(pco)
  dim(pco$li)
  
  ggplot(data=pco$li, aes(x=A1, y=-A2, col=pops[pops%in%humb]))+
    geom_point(aes(pch=type[pops%in%humb]), size=4, alpha=0.8)+
    #stat_chull(aes(fill=pops[pops%in%humb]), geom="polygon",lwd=0,alpha=0.2)+
    scale_color_manual(breaks=levels(pops[pops%in%humb]), values=c(sample(inferno(9),9),"gray20"))+
    #theme_dark()+
    xlab("PCO1 (19.2%)")+
    ylab("PCO2 (9.6%)")+
    labs(col="Population")+
    ggtitle("GBAS Panel (N=109)")+
    theme(text=element_text(size=16), axis.text=element_blank(), axis.ticks=element_blank(),
          panel.grid=element_blank(), panel.background=element_rect(fill='gray80'))+
    guides(fill="none", shape="none")+
    coord_fixed()
  
}



#### PREPARE COORDINATE DATA ####
# fix samples names
samples <- read.csv("data/GBAS_Sample_Information.csv")

# grab mean coordinates for each sample
xs <- aggregate(USAContigConic_x~IndID, samples, FUN=mean)
ys <- aggregate(USAContigConic_y~IndID, samples, FUN=mean)
indcoords <- cbind(xs, ys$USAContigConic_y)
names(indcoords)[2:3] <- c("x","y")
plot(indcoords$x, indcoords$y)#check

# add plot ID
indcoords$PlotID <- inds$PlotID[match(indcoords$IndID, inds$IndID)]
indcoords <- indcoords[-which(indcoords$IndID%in%c("-","FAILED","Unk_F20?","Unk_F27?","Unk_F29?","Unk_M01?","Unk_M08?","Unk_M21?","Unk_M22?")),]
indcoords[indcoords$PlotID %in% c(NA,""),]
indcoords$PlotID[indcoords$PlotID %in% c(NA,"")] <- c("Border_F01.rk","Trinidad_M01.tr","NCali_M05.tr","SOreg_M19.tr",
                                               "NDunes_F01.tr","NDunes_F03.tr","Unk_F30",NA,"NDunes_M03.tr")#unk_M10 is SOreg_M17.tr
length(indcoords$PlotID)==length(unique(indcoords$PlotID))
dim(indcoords)#143 samples with coordinate info

## check samples with missing spatial info
nocoords <- indNames(snpgenind)[which(!indNames(snpgenind) %in% indcoords$PlotID[!is.na(indcoords$x)])]
View(samples[samples$PlotID%in%nocoords, c("SampleID","IndID","PlotID","Easting","Northing","Longitude","Latitude")])
sdn_m16 <- data.frame(IndID="Unk_M16",x=-2237589, y=856711.8, PlotID="SDunes_M16.sc")
indcoords <- rbind(indcoords, sdn_m16)
#indcoords[indcoords$PlotID=="SOreg_M20.tr",c("x","y")] <- c()#missing!

indNames(snpgenind)[!indNames(snpgenind) %in% indcoords$PlotID]
#Karuk hair samples: NCali_M11.sc, NCali_M12.sc
#missing coordinates: NDunes_F10.sc, NDunes_F11.sc, NDunes_M06.sc, NDunes_F12.sc, MACA_HU_F05, BMT_M05, BMT_M09, MAAM_AK.sc

## subset to only samples with spatial info
indcoordsub <- indcoords[which(indcoords$PlotID %in% indNames(snpgenind)),]
dim(indcoordsub)#140 martens (including 1 americana)
snpgenind_sub <- snpgenind[match(indcoordsub$PlotID, indNames(snpgenind))]
snpgenind_sub
sum(!indNames(snpgenind_sub)==indcoordsub$PlotID)#check

## subset to COMA
popsub <- as.character(pop(snpgenind_sub))
humb <- c("Ashland","Trinidad","S Oregon","N California","Border",
  "North Coast", "S Coast OR","S Dunes","N Dunes")
comapops <- which(popsub%in%humb)
coma_coords <- indcoordsub[comapops,]; dim(coma_coords)#100 martens
coma_genind <- snpgenind_sub[match(coma_coords$PlotID, indNames(snpgenind_sub)),]
sum(!indNames(coma_genind)==coma_coords$PlotID)#check



#### PLOT PCOA AXES VS. GEOGRAPHIC SPACE RANGEWIDE ####
# rerun PCoA (full dataset)
pco <- dudi.pco(dist(tab(snpgenind)), nf=2, scannf=FALSE)

# set up df
indcoordsub$PCO1 <- pco$li[match(indcoordsub$PlotID, indNames(snpgenind)), 1]
indcoordsub$PCO2 <- pco$li[match(indcoordsub$PlotID, indNames(snpgenind)), 2]

# function for making interpolation grid
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

# interpolate across grid
img <- interp(indcoordsub$x, indcoordsub$y, indcoordsub$PCO1, duplicate="median", xo=makeGrid(indcoordsub$x), yo=makeGrid(indcoordsub$y))

# load map for plotting
states <- map_data("state",)
statesp <- SpatialPointsDataFrame(states[,c("long","lat")], states, proj4string=CRS("EPSG:4326"))
statesp <- spTransform(statesp, CRS("ESRI:102005"))
states$long <- coordinates(statesp)[,1]
states$lat <- coordinates(statesp)[,2]

# plot PCO1 over map with contours for interpolation surface
ggplot()+
  geom_map(data=states, map=states, aes(map_id=region), fill=NA, col='black')+
  geom_contour_filled(aes(x=rep(img$x, nrow(img$z)), 
                          y=rep(img$y, each=ncol(img$z)),
                          z=as.numeric(img$z)), alpha=0.5, show.legend=FALSE, bins=8)+
  geom_point(data=indcoordsub, aes(x, y, col=PCO1), alpha=0.8, size=3)+
  theme_bw()+
  labs(col="PCO1")+
  coord_fixed()+
  scale_color_viridis(direction=-1)+
  scale_fill_viridis(discrete=T,direction=-1)+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), 
        panel.grid=element_blank(), text=element_text(size=14),
        legend.position="inside", legend.position.inside=c(0.9,0.8))+
  coord_sf(crs=CRS("ESRI:102005"),xlim=c(-2350546.916,-840711.114), ylim=c(200000,1352610.596))

# plot PCO2
# img <- akima::interp(indcoordsub$x, indcoordsub$y, indcoords$PCO2, duplicate="median", xo=makeGrid(indcoordsub$x), yo=makeGrid(indcoordsub$y))
# ggplot()+
#   geom_map(data=states, map=states, aes(map_id=region), fill=NA, col='black')+
#   geom_contour_filled(aes(x=rep(img$x, nrow(img$z)), 
#                           y=rep(img$y, each=ncol(img$z)),
#                           z=as.numeric(img$z)), alpha=0.5, show.legend=FALSE, bins=8)+
#   geom_point(data=indcoordsub, aes(x, y, col=PCO2), alpha=0.8, size=3)+
#   theme_bw()+
#   labs(col="PCO2")+
#   coord_fixed()+
#   scale_color_viridis()+
#   scale_fill_viridis(discrete=T)+
#   theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), 
#         panel.grid=element_blank(), text=element_text(size=14),
#         legend.position="inside", legend.position.inside=c(0.9,0.8))+
#   coord_sf(crs=CRS("ESRI:102005"),xlim=c(-2350546.916,-840711.114), ylim=c(200000,1352610.596))



#### PLOT PCO2 VS. GEOGRAPHIC SPACE HUMBOLDT ####
# rerun pco
pco <- dudi.pco(dist(tab(snpgenind[pop(snpgenind)%in%humb,])), nf=2, scannf=FALSE)
coma_coords$PCO1 <- pco$li[match(coma_coords$PlotID, rownames(pco$tab)), 1]
coma_coords$PCO2 <- pco$li[match(coma_coords$PlotID, rownames(pco$tab)), 2]

# convert coords to UTM
sp <- SpatialPointsDataFrame(coma_coords[,c("x","y")], data=coma_coords, proj4string=CRS("ESRI:102005"))
utm <- spTransform(sp, CRS("EPSG:32610"))
utm$x <- utm@coords[,1]
utm$y <- utm@coords[,2]

# convert basemap to UTM
statesp <- SpatialPointsDataFrame(states[,c("long","lat")], states, proj4string=CRS("EPSG:4326"))
#statesp <- spTransform(statesp, CRS("EPSG:32610"))


# plot PCO1
img <- interp(utm$x, utm$y, utm$PCO2, duplicate="median", xo=makeGrid(utm$x), yo=makeGrid(utm$y))
ggplot()+
  geom_map(data=states, map=states, aes(map_id=region), fill=NA, col='black')+
  geom_contour_filled(aes(x=rep(img$x, nrow(img$z)), 
                          y=rep(img$y, each=ncol(img$z)),
                          z=as.numeric(img$z)), alpha=0.5, show.legend=FALSE, bins=8)+
  geom_point(data=utm@data, aes(x, y, col=PCO2), alpha=0.8, size=3)+
  theme_bw()+
  labs(col="PCO2")+
  #coord_fixed()+
  scale_color_viridis()+
  scale_fill_viridis(discrete=T)+
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), 
        text=element_text(size=14),
        legend.position="inside", legend.position.inside=c(0.8,0.8))+
  coord_sf(crs=CRS("EPSG: 32610"),xlim=c(367518,492244), ylim=c(4524138,4934554))#extent of QGIS MAP


#####################################
#### CALCULATE FST & RELATEDNESS ####
#####################################
# pairwise Fst
hierfstat::pairwise.neifst(snpgenind_sub)
hierfstat::boot.ppfst(snpgenind_sub, nboot=999)

# calculate relatedness
betas <- as.data.frame(hierfstat::beta.dosage(snpgenind_sub))
# factor pop for proper coloring...
betas$Pop <- unlist(lapply(colnames(betas), function(x) strsplit(x, split="_")[[1]][1]))
betas$Pop <- factor(betas$Pop, levels=c("americana","NCoast","NDunes","SDunes","SCoast","SOreg","Applegate",
                                        "Border","NCali","MapleCrk","ORcasc","Lssn","ORbmt","WAcasc","WY","CO","MACA"))
rownames(betas)[is.na(betas$Pop)]
betas$Pop[is.na(betas$Pop)] <- "ORcasc"
# plot for key samples
ggplot(betas[!betas$Pop%in%c("BMT","CO","MAAM","WACasc","WY","Lssn","ORCasc"),])+
  geom_point(aes(x="Applegate_M01.sc", y=Applegate_M01.sc, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  geom_point(aes(x="NCoast_F01.sc", y=NCoast_F01.sc, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  geom_point(aes(x="NCoast_F02.sc", y=NCoast_F02.sc, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  geom_point(aes(x="NCoast_F03.sc", y=NCoast_F03.sc, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  geom_point(aes(x="NCoast_F04.sc", y=NCoast_F04.sc, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  geom_point(aes(x="SOreg_M10.tr", y=SOreg_M10.tr, col=Pop), size=2, alpha=0.5, position=position_jitter(height=0, width=0.3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  xlab(element_blank())+ylab("Pairwise Relatedness")+
  labs(col="Population")+
  scale_color_manual(values=c("white",rev(pal)))



#######################################
#### CONSTRUCT TREES FOR GBAS DATA ####
#######################################
# remove samples with <80 SNPs
remove <- indNames(snpgenind)[apply(snpgenind@tab, 1, function(x) sum(is.na(x)))>((94-75)*2)]; remove
snpgenind_sub <- snpgenind[-which(indNames(snpgenind)%in%remove),]
indNames(snpgenind_sub)[which(indNames(snpgenind_sub)=="MAAM_WA.tr")] <- "americana_WA.tr"
indNames(snpgenind_sub)[which(indNames(snpgenind_sub)=="MAAM_AK.sc")] <-"americana_AK.sc"
gipops <- unlist(lapply(indNames(snpgenind_sub), function(X) strsplit(X,"_")[[1]][1]))
gipops[which(indNames(snpgenind_sub)=="MACA_HU_F05")] <- "Unk"
gipops <- factor(gipops, levels=c("americana","WY","CO","WACasc","BMT","Lssn","ORCasc","Trinidad","NCali",
                                  "Border","Ashland","SOreg","SCoast","SDunes","NDunes","NCoast","Unk"))


# calculate prevosti distance- fraction of dissimilar sites
prev <- poppr::prevosti.dist(snpgenind_sub@tab)
#bitd <- poppr::bitwise.dist(gl)#this is basically euclidean distance
gl <- as.genlight(snpgenind_sub@tab)
euc <- dist(gl)

## Function for color-coded plotting
plotTree <- function(tree, genind=snpgenind_sub, initpops=gipops, title=NULL){
  # ladderize
  tree <- tree %>% ladderize()
  # root tree
  roottree <- root(tree, outgroup=c("americana_WA.tr","americana_AK.sc"), resolve.root=TRUE)
  # set up pops in order for plotting
  pops <- initpops[match(roottree$tip.label, indNames(genind))]
  # plot
  myPal <- c('black',rev(rainbow(15,v=0.5)))
  par(mar=c(2,1,1,1))
  plot(roottree, show.tip=FALSE)
  title(title)
  tiplabels(roottree$tip.label, col=myPal[pops], bg=NULL, frame="none",cex=0.5, fg="transparent",adj=c(0,0.05))
  #axisPhylo()
}#plotTree


#### UPGMA TREE ####
marten_UPGMA <- upgma(euc)
plotTree(marten_UPGMA)


#### NEIGHBOR-JOINING TREE ####
marten_NJ <- njs(euc)
marten_NJroot <- root(marten_NJ, outgroup=c("americana_WA.tr","americana_AK.sc"))
# plot
plotTree(marten_NJ)
# assess certainty
bstrees <- boot.phylo(marten_NJroot, snpgenind_sub@tab, njs, trees=TRUE, B=100)#make 100 bootstrap trees
clad <- prop.clades(marten_NJroot, bstrees$trees, rooted=FALSE)#count # trees each bipartition is present in
hist(clad)
# add arbitrary node labels (otherwise won't open in FigTree)
marten_NJroot <- ape::makeNodeLabel(marten_NJroot); marten_NJroot


#### MAXIMUM-LIKELIHOOD TREE ####
# perform model selection for substition model
mt <- modelTest(phydat, model=c("JC","F81","GTR","SYM","TVM","HKY"))
# check fit of NJ tree
dist<- dist.ml(phydat)
fit <- pml(marten_NJroot, phydat)
print(fit)
# run ML tree (rates and edges can't be optimized at same time)
fitML <- optim.pml(fit, optNni=TRUE)#, model="JC")
                   #optNni=TRUE,#optimize tree topology
                   #optEdge=TRUE),#optimize edge lengths
                   #optRate=TRUE,#optimize overall rate
                   #optBf=TRUE,#optimize base freqs
                   #optQ=TRUE,#optimize rate of possible substitutions
                   #optGamma=TRUE,#use gamma dist to model rate of substitution across sites
                   #ASC=TRUE)#use ascertainment bias correction models
# check fits
logLik(fit)
logLik(fitML)
anova(fit, fitML)
AIC(fit)
AIC(fitML)
# extract new tree
MLroot <- root(fitML$tree, outgroup=c("americana_AK.sc","americana_WA.tr"))
MLroot <- ladderize(MLroot)
# calculate bootstrap support
bs <- bootstrap.pml(fitML, bs=100, ASC=TRUE, optNni=TRUE, optEdge=TRUE, #optRate=TRUE,optBf=TRUE,optQ=TRUE,optGamma=TRUE,
                    multicore=FALSE, control = pml.control(trace=0))
# plot bootstraps
plotBS(MLroot, bs, p=0, type="p", cex=0.5, bs.col="red")
axisPhylo()

plotTree(MLroot)
axisPhylo()
#add.scale.bar(x=-1,length=10)

#poppr::aboot(snptab, dist=poppr::prevosti.dist, sample=100, tree="nj", cutoff=50)
#drawSupportonEdges
plot(marten_NJ, cex=0.6, col=rainbow(17)[pops])
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

groups <- split(marten_NJ$tip.label, gsub("_\\w+", "", marten_NJ$tip.label))
#groups <- factor(groups, levels=levels(pops))
martenNJ_groups <- groupOTU(marten_NJ, groups)

# ggtree
ggtree(martenNJ_groups, aes(col=group), right=FALSE, layout="dendrogram", size=1.5)+#layout="fan", )+
  geom_treescale(y=-5, offset=-3, fontsize=5)+
  scale_color_manual(values=c(rainbow(14,v=0.9)[c(1:9,9,10,11,11:14)],'grey'),#[c(8,14,10,6,1,2,12,9,9,11,4,3,5,7,11,13)]))+
                     breaks=c("NCoast","NDunes","SDunes","SCoast","SOreg","NCali","Trinidad","Ashland",
                              "ORcasc","ORcasc-2","Lssn","WAcasc","Rainier","ORbmt","WY","CO"),
                     labels=c("N Coast","N Dunes","S Dunes","S Coast","S Oregon","N California","Trinidad","Ashland",
                              "OR Cascades-1", "OR Cascades-2","Lassen"," WA Cascades","Mt Rainier","Blue Mtns","Wyoming","Colorado"))+
  geom_tiplab(size=3)+
  #scale_x_continuous(expand=c(0.001,0))+
  theme_tree()+
  theme(plot.margin=unit(c(5,5,20,5),"mm"), text=element_text(size=18))+
  labs(col="Population")+
  xlab("Prevosti distance")



###############################################
#### CALCULATE PID & PIDsib PER POPULATION ####
###############################################
# function to calculate PIDsib based on MAFs
# Waits et al. 2001 function
PIDsib <- function(maf){
# mafs = list of minor allele frequencies, 1 value per bilallelic locus
  pi <- maf
  pj <- 1-maf
  sumsq <- pi^2 + pj^2 
  sumqt <- pi^4 + pj^4
  PIDsib_locus <- 0.25 + (0.5*sumsq) + (0.5*(sumsq^2)) - (0.25*sumqt)
  PIDsib <- prod(PIDsib_locus)
  return(PIDsib)
}#PIDsib


# function to calculate multilocus PID
PID <- function(maf){
  # maf = minor allele frequency (per locus)
  pi <- maf # probability of allele 1 (ref)
  pj <- 1-maf # probability of allele 2 (alt)
  PID_locus <- pi^4 + pj^4 + (2*pi*pj)^2 # combine per locus calc
  PID <- prod(PID_locus) # take product for multilocus
  return(PID)
}#PID


# function to calculate PID based on MAFs
# based on Paetkau and Strobeck 1994, in Waits et al. 2001
PIDunbiased <- function(maf, n){
  # maf = minor allele frequency (per locus)
  pi <- maf
  pj <- 1-maf
  a2 <- pi^2 + pj^2
  a3 <- pi^3 + pj^3
  a4 <- pi^4 + pj^4
  PIDunbiased_locus <- ((n^3)*(2*a2^2-a4)-2*(n^2)*(a3+2*a2)+n*(9*a2+2)-6) /
    ((n-1)*(n-2)*(n-3))
  PIDunbiased <- prod(PIDunbiased_locus)
  return(PIDunbiased)
}#PIDunbiased


# wrapper function to calculate PID, PIDsib, PIDunbiased per row in genpop
AssessSNPpanelByPop <- function(genind, threshold=0){
  # genind : genind object with populations specified
  # threshold: MAF threshold for counting & including SNPs
  # set up output DF
  outDF <- data.frame()
  # run calcs for each population
  for (pop in levels(pop(genind))){
    sub <- genind[pop(genind)==pop,] # subset genind
    nsamp <- nrow(sub@tab) # calc number samples
    mafs <- unlist(minorAllele(sub)) # calculate minor allele freqs
    mafs <- mafs[mafs>threshold] # remove SNPs below threshold
    nsnps <- length(mafs) # count SNPs remaining
    PIDsib <- PIDsib(mafs) # run PIDsib
    PID <- PID(mafs) # run PID
    PIDunbiased <- PIDunbiased(mafs, nsamp) # run PIDunbiased
    row <- c(pop, nsamp, nsnps, PIDsib, PID, PIDunbiased)# make row for output
    outDF <- rbind(outDF, row) # add row to output DF
  }#pop
  # clean up output
  names(outDF) <- c("Population","N_Samples","N_SNPs","PIDsib","PID","PIDunbiased")
  return(outDF)
  
}#AssessSNPpanelByPop


# run calcs for each population and for montane vs. humboldt
AssessSNPpanelByPop(snpgenind)

temp <- snpgenind
pops <- as.character(pop(temp))
pops[pops%in%c("North Coast","N Dunes","S Dunes","S Coast","Ashland","S Oregon","N California","Border","Trinidad","Unknown")] <- "Humboldt"
pops[pops%in%c("OR Cascades","Lassen","WA Cascades","Wyoming","Colorado","Blue Mountains")] <- "Mountain"
pop(temp) <- pops
AssessSNPpanelByPop(temp)

