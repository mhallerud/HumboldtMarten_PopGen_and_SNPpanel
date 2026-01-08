########## COMPARISON FOR POPULATION GENETIC INFERENCE BETWEEN RADSEQ & GBAS SNP PANEL ############
#### AUTHOR: Maggie Hallerud
#### PAPER: Genomics and noninvasive genetics quantify metapopulation structure and genetic viability of an elusive and threatened small carnivore

#setwd('/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/RADseq_Plate3_Results/')
library(vcfR)
library(adegenet)
library(dartR)
library(radiator)
library(ASRgenomics)
library(SeqArray)
library(stringr)
library(SNPRelate)
library(hierfstat)


#### LOAD AND FILTER RADSEQ DATASET ####
# input RADseq LD-pruned dataset
radvcf <- vcfR::read.vcfR('LDprunedRAD.vcf')
fixes <- as.data.frame(radvcf@fix)
radvcf@fix[,1] <- unlist(lapply(fixes$ID, function(x) strsplit(x, ":")[[1]][1]))

# convert to genind & replace names & add pop
genind <- vcfR::vcfR2genind(radvcf)
match(indNames(genind), popmap$MergeID)
pop(genind) <- unlist(lapply(indNames(genind), function(x) strsplit(x, "_")[[1]][1]))



#### LOAD SNP PANEL DATASET ####
## load CSV formatted SNPs
samples <- read.csv("GBAS_sample_identification.csv")

# convert to plink-formatted genotypes
{
  inds <- samples[which(samples$Use==1),]

  genos <- inds[,c(13:115)]
  genos[genos=="0"] <- NA
  genos[genos=="NA"] <- NA
  genos[genos=="-"] <- NA
  genos[genos=="FALSE"] <- NA
  sort(unique(unlist(genos)))
  inds[,13:115] <- genos
  
  inds <- as.data.frame(inds)
  
  ## convert to A1HOM / A2HOM / HET
  dim(inds)
  plink <- inds
  # for (j in 11:109){
  #   # grab alleles
  #   calls <- sort(unique(inds[,j]))
  #   calls <- calls[!is.na(calls)]
  #   i <- 0
  #   for (c in calls){
  #     i=i+1
  #     # identify A1HOM, A2HOM, and HET genotypes
  #     alleles <- strsplit(c,"")[[1]]
  #     if((alleles[1]==alleles[2]) & i==1){
  #       a1hom <- c
  #     }else{
  #       if((alleles[1]==alleles[2]) & i>1){
  #         a2hom <- c
  #       }else{
  #         het <- c
  #       }#ifelse
  #     }#ifelse
  #     # rename genotypes as A1HOM/A2HOM/HET
  #     plink[which(inds[,j]==a1hom), j] <- "A1HOM"
  #     plink[which(inds[,j]==a2hom), j] <- "A2HOM"
  #     plink[which(inds[,j]==het), j] <- "HET"
  #   }#c
  # }#j
  
  # now convert to plink format
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
  genos[genos=="GG* (1 rep)"] <- "0 0"
  genos[is.na(genos)] <- "0 0"
  plink <- inds
  plink[,13:115] <- genos
  
  names(plink) <- gsub("\\.", "_",names(plink))
  plink$Pop[is.na(plink$Pop)]#check
  #View(plink)#check
  #write.csv(plink, "individualIDs_17jun2025_plink.csv")
}

# read in PLINK format data and convert to genind for processing
#plink <- read.table("ALL_Genotypes_Plink_Format.raw", sep="\t")
length(plink$PlotID)#150 individuals
length(unique(plink$PlotID))
plink$PlotID[plink$PlotID=="MAAM"] <- c("americana_WA.tr","americana_AK.sc")
snpgenind <- adegenet::df2genind(plink[,13:115], sep=" ", NA.char="0 0", ind.names=plink$PlotID, pop=plink$Pop, ncode=3)
snpgenind#150 martens x 103 loci (including SRY! + microhaps)
#col51 is nonspecific MACA_67356.0:111

# convert and write to VCF
#radiator::genomic_converter(snpgenind, filename="GBAS_individuals.vcf", output=c("vcf","plink"))
#snpvcf <- read.vcfR("GBAS_individuals.vcf")
#snpvcf@fix[,3] <- snpvcf@fix[,2] #CHROM=POS (POS WAS LABELED BY AMPLICON)
#gt <- as.data.frame(vcf@gt)

# extract populations from individual IDs
indNames(snpgenind)[indNames(snpgenind)=="MACA_HU_F05"] <- "UNK_MACAHUF05"
pop(snpgenind) <- unlist(lapply(indNames(snpgenind), function(x) strsplit(x, "_")[[1]][1]))

# correct indNames so that they match the RAD data
indNames(snpgenind) <- unlist(lapply(indNames(snpgenind), function(x) strsplit(x,"\\.")[[1]][1]))
# fix different names
indNames(snpgenind)[indNames(snpgenind)=="NCali_M04"] <- "NCali_Mrdkl1_2"
indNames(snpgenind)[indNames(snpgenind)=="NCali_M06"] <- "NCali_Mrdkl2"
indNames(snpgenind)[indNames(snpgenind)=="SOreg_M18"] <- "SOreg_M22" #M18 and M22 are same animal
indNames(snpgenind)[indNames(snpgenind)=="Border_F01"] <- "Border_F01rk"


#### SUBSET DATASETS FOR MATCHING INDIVIDUALS ####
sum(indNames(genind)%in%indNames(snpgenind))#59
indNames(genind)[!indNames(genind)%in%indNames(snpgenind)]#check excluded

radsamps <- which(indNames(genind) %in% indNames(snpgenind))
genind_sub <- genind[radsamps,]

snpsamps <- match(indNames(genind_sub),indNames(snpgenind))
snpgenind_sub <- snpgenind[snpsamps,]

#check...
sum(indNames(genind_sub)!=indNames(snpgenind_sub))



#### TISSUE SAMPLES GENOTYPED BY GBAS BUT NOT RADSEQ ####
# extract sample types
sampletype <- unlist(lapply(plink$PlotID, function(x) strsplit(x,"\\.")[[1]][2]))
plink$PlotID[is.na(sampletype)]
sampletype[is.na(sampletype)] <- "rk" #fix NAs
sampletype_sub <- sampletype[snpsamps]

# check which roadkill and trapping samples were not included in RADseq
unique(sampletype)
indNames(snpgenind)[which(!indNames(snpgenind)%in%indNames(genind) & sampletype%in%c("rk","tr"))]
# note: one of these is american marten (G40404 / MAAM_WA)



########################
#### ROUND 1: PCO ######
########################

# RAD data 
radtab <- tab(genind_sub, freq=TRUE, NA.method="asis")
radpco <- dudi.pco(dist(radtab), scannf=FALSE, nf=54)
s.class(radpco$li, fac=pop(genind_sub))
screeplot(radpco)
summary(radpco)

# SNP data
snptab <- tab(snpgenind_sub, freq=TRUE, NA.method="asis")
snppco <- dudi.pco(dist(snptab), scannf=FALSE, nf=54)
s.class(snppco$li, fac=pop(snpgenind_sub))
screeplot(radpco)

# pop names for plotting
pops <- as.character(pop(snpgenind_sub)); unique(pops)
pops[pops=="NDunes"] <- "N Dunes"
pops[pops=="SDunes"] <- "S Dunes"
pops[pops=="SOreg"] <- "S Oregon"
pops[pops=="NCali"] <- "N California"
pops[pops=="MapleCrk"] <- "Maple Creek"
pops[pops=="Lssn"] <- "Lassen"
pops[pops=="ORcasc"] <- "OR Cascades"
pops[pops=="WAcasc"] <- "WA Cascades"
pops[pops=="ORbmt"] <- "Blue Mtns"
pops[pops=="CO"] <- "N Colorado"
pops[pops=="WY"] <- "Wyoming"
unique(pops)

# factor pops including levels that are present only in GBAS s.t. colors are consistent
pops <- factor(pops, levels=c("N Coast","N Dunes","S Dunes","S Oregon","S Coast",
                              "Border","Maple Creek","N California","Applegate",
                              "OR Cascades","Lassen","WA Cascades","Blue Mtns",
                              "Wyoming","Colorado","Unk"))

cols <- c("#000000",#black-UNK
  "#db6d00",#orange-CO
  "gray20",#"#252525",#darkgray-WY
  "gray60",#676767",#lightgray-ORbmt
  "gray90",#"#ffffff",#white-WAcasc
  "mediumpurple4",#171723",#navy-Lssn
  "darkcyan",#"#004949",#dark teal-ORcasc
  "aquamarine1",#009999",#light teal-Applegate
  "#22cf22",#bright green-NCali
  "darkviolet",#"#490092",#purple-MapleCrk
  "steelblue1",#"#006ddb",#sky blue-Border
  "#b66dff",#light purple-SCoast
  "#ff6db6",#bright pink-SOreg
  "#920000",#dark red-SDunes
  "#ffdf4d",#yellow-NDunes
  "#8f4e00")#brown-NCoast

PC1 <- ggplot()+
  geom_point(aes(x=radpco$li$A1, y=radpco$li$A2, col=pops), size=3, alpha=0.7)+
  theme_classic()+
  xlab("PCO1 (18.0%)")+
  ylab("PCO2 (7.0%)")+
  ggtitle("ddRADseq")+
  scale_color_manual(values=rev(cols)[levels(pops)%in%pops], breaks=levels(pops)[levels(pops)%in%pops])+
  theme(legend.position="none",text=element_text(size=16), 
        axis.text=element_blank(), axis.ticks=element_blank())
PC2 <- ggplot()+
  geom_point(aes(x=snppco$li$A1, y=snppco$li$A2, col=pops), size=3, alpha=0.7)+
  theme_classic()+
  xlab("PCO1 (21.1%)")+
  ylab("PCO2 (12.1%)")+
  ggtitle("GBAS")+
  scale_color_manual(values=rev(cols)[levels(pops)%in%pops], breaks=levels(pops)[levels(pops)%in%pops])+
  theme(legend.position="none",text=element_text(size=16), 
        axis.text=element_blank(), axis.ticks=element_blank())

ggpubr::ggarrange(PC1, PC2, widths=c(1,1))



################################################
#### ROUND 2: INDIVIDUAL SNP HETEROZYGOSITY ####
################################################
# calculate heterozygosity for each
radhet <- apply(radtab, 1, function(X) sum(X=="0.5",na.rm=T)/sum(!is.na(X)))
hist(radhet)

snphet <- apply(snptab, 1, function(X) sum(X=="0.5",na.rm=T)/sum(!is.na(X)))
hist(snphet)

# check correlations
cor(snphet, radhet, method="pearson")#-0.10
cor(snphet, radhet, method="spearman")#-0.24

plot(radhet, snphet, col=cols[pops], pch=16)


## Try again- filter RAD data for MAF>0.05 and 1 SNP per locus
radmafs <- maf(radvcf[,c(1,radsamps+1)], element=2)
radmafs[radmafs>50] <- 100-radmafs[radmafs>50]
radmafs <- as.data.frame(radmafs)
radmafs$ID <- rownames(radmafs)
head(radmafs)

fixes$CHROM <- unlist(lapply(fixes$ID, function(X) strsplit(X,":")[[1]][1]))
singleSNPs <- unlist(lapply(unique(fixes$CHROM), FUN=function(x){
  pos <- fixes$ID[fixes$CHROM==x]
  snps <- radmafs[radmafs$ID%in%pos,]
  snp <- sample(snps$ID[snps$Frequency==max(snps$Frequency)], 1)
  return(snp)
}))#lapply
length(singleSNPs)
length(unique(fixes$CHROM))

mafs05 <- radmafs$ID[radmafs$Frequency>0.05]

radgt <- as.data.frame(extract.gt(radvcf)[,c(radsamps)])
radgt_maf05 <- radgt[which(rownames(radgt)%in%singleSNPs & rownames(radgt)%in%mafs05),]#7002
radmaf5_het <- apply(radgt_maf05, 2, function(X) sum(X=="0/1",na.rm=T)/sum(!is.na(X)))


cor(radmaf5_het, snphet, method="pearson")#0.04
cor(radmaf5_het, snphet, method="spearman")#0.001
plot(radmaf5_het, snphet)


## Final attempt: try LD-pruning in addition to the above 
radtab <- tab(genind_sub, freq=F, NA.method="asis")
radtab_filt <- ASRgenomics::qc.filtering(M=radtab, maf=0.01)
rad_LDpruned <- ASRgenomics::snp.pruning(M=radtab_filt$M.clean, method="correlation", pruning.thr=0.60)
dim(rad_LDpruned$Mpruned)#8308

radprune_het <- apply(rad_LDpruned$Mpruned, 1, function(X) sum(X=="1",na.rm=T)/sum(!is.na(X)))

cor(radprune_het, snphet, method="pearson")#0.11
cor(radprune_het, snphet, method="spearman")#0.11
plot(radprune_het, snphet)

## Presenting unfiltered results since it doesn't seem to matter...
ggplot()+
  geom_abline(aes(slope=1, intercept=0))+
  geom_point(aes(x=radhet, y=snphet, col=pops), alpha=0.8, size=3)+
  theme_classic()+
  scale_color_manual(values=rev(cols), breaks=levels(pops))+
  ggtitle(expression(H[SNP]))+
  xlab("ddRADseq")+
  ylab("GBAS")+
  theme(legend.position="none", text=element_text(size=18))
## save as 450 x 400



#########################
#### ROUND 3: KINSHP ####
#########################
# first we'll convert everything to SNPGDS format...
SNPRelate::snpgdsCreateGeno(gds.fn="RAD_LDpruned.gds", 
                            genmat=tab(genind_sub), sample.id=indNames(genind_sub), snpfirstdim=FALSE)
rad_gds <- SNPRelate::snpgdsOpen("RAD_LDpruned.gds")

SNPRelate::snpgdsCreateGeno("SNPpanel_RADsubset_95plex.gds", 
                            genmat=tab(snpgenind_sub), sample.id=indNames(snpgenind_sub), snpfirstdim=FALSE)
snp_gds <- SNPRelate::snpgdsOpen("SNPpanel_RADsubset_95plex.gds")

# RAD
rad_betas <- SNPRelate::snpgdsIndivBeta(rad_gds, sample.id=read.gdsn(index.gdsn(rad_gds, "sample.id")), inbreeding=TRUE, with.id=TRUE, )#diaagonal: individual variance est.
rad_betasc <- SNPRelate::snpgdsIndivBetaRel(rad_betas, min(rad_betas$beta))#correct beta estimates
rownames(rad_betas$beta) <- rad_betas$sample.id
colnames(rad_betas$beta) <- rad_betas$sample.id
corrplot(rad_betas$beta)
rad_betas_long <- reshape2::melt(rad_betas$beta)
rad_betasc_long <- reshape2::melt(rad_betasc$beta)
#rad_betasc_long <- rad_betasc_long[-which(rad_betasc_long$Var1==rad_betasc_long$Var2),]#remove self-comparisons


# SNP 95-plex raw
snp_betas <- SNPRelate::snpgdsIndivBeta(snp_gds, sample.id=read.gdsn(index.gdsn(snp_gds, "sample.id")), inbreeding=FALSE, with.id=TRUE)#diaagonal: individual variance est.
snp_betasc <- SNPRelate::snpgdsIndivBetaRel(snp_betas, min(snp_betas$beta))
rownames(snp_betas$beta) <- indNames(snpgenind_sub)
colnames(snp_betas$beta) <- indNames(snpgenind_sub)
corrplot(snp_betas$beta)
snp_betas_long <- reshape2::melt(snp_betasc$beta)
#snp_betas_long <- snp_betas_long[-which(snp_betas_long$Var1==snp_betas_long$Var2),]#remove self-comparisons

# check correlations
cor(rad_betasc_long$value, snp_betas_long$value, method="pearson")#0.59
cor(rad_betasc_long$value, snp_betas_long$value, method="spearman")#0.50
plot(rad_betasc_long$value, snp_betas_long$value, col=cols[pops], pch=16)

# check correlations within Humboldt martens only
rad_betas_long$Pop1 <- unlist(lapply(rad_betas_long$Var1, function(X) strsplit(as.character(X),"_")[[1]][1]))
rad_betas_long$Pop2 <- unlist(lapply(rad_betas_long$Var2, function(X) strsplit(as.character(X),"_")[[1]][1]))
humboldt <- c("SDunes","NDunes","NCali","Trinidad","SOreg")
betas_humboldt <- which(rad_betas_long$Pop1 %in% humboldt & rad_betas_long$Pop2 %in% humboldt)
cor(rad_betasc_long$value[betas_humboldt], snp_betas_long$value[betas_humboldt], method="pearson")#0.69
cor(rad_betasc_long$value[betas_humboldt], snp_betas_long$value[betas_humboldt], method="spearman")#0.62
cor(rad_betasc_long$value[betas_humboldt&which(rad_betas_long$Pop1==rad_betas_long$Pop2)], 
    snp_betas_long$value[betas_humboldt&which(rad_betas_long$Pop1==rad_betas_long$Pop2)], method="pearson")#0.69


# check correlations within each pop
betasWithin <- which(rad_betas_long$Pop1==rad_betas_long$Pop2)
cor(rad_betasc_long$value[betasWithin], snp_betas_long$value[betasWithin], method="pearson")#0.49
cor(rad_betasc_long$value[betasWithin], snp_betas_long$value[betasWithin], method="spearman")#0.37
plot(rad_betasc_long$value[betasWithin], snp_betas_long$value[betasWithin])

cor(rad_betasc_long$value[c(betasWithin,betas_humboldt)], snp_betas_long$value[c(betasWithin,betas_humboldt)])


# plot results - color by within-group comparisons
df <- data.frame(rad_betas_long$Pop1, rad_betas_long$Pop2, rad_betas_long$Var1, rad_betas_long$Var2,
                 rad_betasc_long$value, snp_betas_long$value)
names(df) <- c("Pop1","Pop2","ID1","ID2","RADbeta","SNPbeta")
df$Pop1[df$Pop1=="Border"] <- "NCali" #bin Border w/ NCali for simplicity & because groups with NCali genetically
df$Pop2[df$Pop2=="Border"] <- "NCali"
df$Pop1[df$Pop1%in%c("ORbmt")] <- "BlueMts"
df$Pop2[df$Pop2%in%c("ORbmt","WAcasc","CO","WY")] <- "E Montane"
df$Pop2[df$Pop2%in%c("ORcasc","Lssn")] <- "W Montane"
df$Pop2 <- factor(df$Pop2, levels=c("NDunes","SDunes","SOreg","NCali","MapleCrk","W Montane","E Montane"))
df$Pop1 <- factor(df$Pop1, levels=c("NCoast","NDunes","SDunes","SOreg","SCoast","Border","MapleCrk",
                                    "NCali","Applegate","ORcasc","Lssn","WAcasc","BlueMts","WY","CO","Unk"))
ggplot(df)+
  geom_abline(aes(slope=1, intercept=0))+
  geom_point(aes(x=RADbeta, y=SNPbeta, col=as.factor(Pop1)), alpha=0.5)+
  facet_wrap(~Pop2)+
  xlab("ddRADseq")+ylab("GBAS")+
  theme_classic()+
  ggtitle(expression(paste("Relatedness (",beta[jj],")")))+
  labs(col="Population")+
  scale_color_manual(values=rev(cols), breaks=levels(df$Pop1))+
  theme(text=element_text(size=16))

df$PopGroup <- paste(df$Pop1, df$Pop2, sep="-")
df$PopGroup[which(df$Pop1!=df$Pop2)] <- df$Pop1[which(df$Pop1!=df$Pop2)]
unique(df$PopGroup)
df$PopGroup <- factor(df$PopGroup, levels=c("NDunes-NDunes","SDunes-SDunes","SOreg-SOreg","NCali-NCali","MapleCrk-MapleCrk",
                                            "ORcasc-ORcasc","Lssn-Lssn","WAcasc-WAcasc","ORbmt-ORbmt","WY-WY","CO-CO",
                                            "NDunes","SDunes","SOreg","NCali","MapleCrk","ORcasc","WAcasc","ORbmt","WY","CO"))
df$Subspecies <- NA
df$Subspecies[df$PopGroup%in%c("NDunes","SDunes","SOreg","NCali","MapleCrk")] <- "Humboldt"
df$Subspecies[df$Subspecies!="Humboldt"] <- "Mountain"
df$Subspecies <- factor(df$Subspecies, levels=c("Humboldt","Mountain"))

df$Pop1 <- factor(df$Pop1, levels=c("NDunes","SDunes","SOreg","NCali","MapleCrk","ORcasc","WAcasc","ORbmt","WY","CO"))
ggplot(df)+
  geom_abline(aes(slope=1, intercept=0))+
  geom_point(aes(x=RADbeta, y=SNPbeta), col=c("gray","black")[df$Subspecies],alpha=0.5, size=1)+
  geom_point(data=df[df$Pop1==df$Pop2,], aes(x=RADbeta, y=SNPbeta, col=Pop1), alpha=0.8, size=2)+
  scale_color_manual(values=cols)+
  theme_classic()+
  theme(legend.position="none", text=element_text(size=16))+
  xlab("ddRADseq")+ylab("GBAS")+ggtitle(expression(beta))



#############################
#### ROUND 4: INBREEDING ####
#############################
# calculate inbreeding coefficients using Weir and Goudet metric
rad_Fbeta <- hierfstat::beta.dosage(genind_sub, Mb=TRUE, inb=TRUE)
rad_Fbeta_inb <- diag(rad_Fbeta$betas)
snp_Fbeta <- hierfstat::beta.dosage(snpgenind_sub, Mb=TRUE, inb=TRUE)
snp_Fbeta_inb <- diag(snp_Fbeta$betas)

# check correlations
cor(rad_Fbeta_inb, snp_Fbeta_inb, method="pearson")#-0.10
cor(rad_Fbeta_inb, snp_Fbeta_inb, method="spearman")#-0.24
plot(rad_Fbeta_inb, snp_Fbeta_inb, pch=16, col=cols[pops])

hom_df <- data.frame(IID=indNames(genind_sub), FID=pop(genind_sub), 
                     RADbeta=rad_Fbeta_inb, SNPbeta=snp_Fbeta_inb)
hom_df$FID <- factor(hom_df$FID, levels=c("NCoast","NDunes","SDunes","SOreg","SCoast","Border","MapleCrk","NCali",
                                          "Applegate","ORcasc","Lssn","WAcasc","ORbmt","WY","CO","Unk"))
ggplot(hom_df)+
  geom_abline(aes(slope=1,intercept=0))+
  geom_point(aes(x=RADbeta, y=SNPbeta, col=FID), size=3, alpha=0.8)+
  theme_bw()+
  ggtitle(expression(F[beta]))+
  coord_flip()+
  xlab("GBAS")+
  ylab("ddRADseq")+
  theme_classic()+
  scale_color_manual(values=rev(cols), breaks=levels(hom_df$FID))+
  labs(color="Population")+
  theme(text=element_text(size=18))
##save as 550 x 400




######################
#### ROUND 5: Fst ####
######################
# calculate Fsts
radhier <- hierfstat::genind2hierfstat(genind_sub)
radfst <- hierfstat::pairwise.neifst(radhier)

snphier <- hierfstat::genind2hierfstat(snpgenind_sub)
snpfst <- hierfstat::pairwise.neifst(snphier)

# corrplot for SNP results
snpfst <- snpfst[c("NDunes","SDunes","SOreg","NCali","Trinidad","ORcasc","Lssn","WAcasc","ORbmt","WY","CO"), 
                 c("NDunes","SDunes","SOreg","NCali","Trinidad","ORcasc","Lssn","WAcasc","ORbmt","WY","CO")]
rownames(snpfst) <- colnames(snpfst) <- c("N Dunes", "S Dunes", "S Oreg", "N Calif", "Trinidad", 
                                          "OR Casc", "Lassen", "WA Casc", "Blue Mtn", "WY", "N CO")
corrplot.mixed(snpfst, lower="number", upper="shade", upper.col=rev(COL2('RdBu',20)), lower.col='black', col.lim=c(0,1), diag='n')

radfst_long <- reshape2::melt(radfst)
snpfst_long <- reshape2::melt(snpfst)
cor(radfst_long$value, snpfst_long$value, method="pearson", use="complete.obs")#0.87
cor(radfst_long$value, snpfst_long$value, method="spearman", use="complete.obs")#0.86

ggplot()+
  geom_point(aes(x=radfst_long$value, y=snpfst_long$value))+
  theme_classic()+
  xlab("ddRADseq")+
  ylab("GBAS")+
  theme(text=element_text(size=16))+
  ggtitle(expression(F[ST]))



######################################################################
#### CALCULATE fORCA BASED ON SNP PANEL RESULTS (ALL INDIVIDUALS) ####
######################################################################

## eqtns found here: https://rosenberglab.stanford.edu/papers/maximin.pdf

## clean up vcf before proceeding...
allpops <- as.character(pop(snpgenind))
#allpops[c(37:39,41)] <- "ORcasc"
#allpops[c(81,86,87)] <- "NCoast"
#allpops[c(82,97,100:101,128,130)] <- "NDunes"
#allpops[c(83:85,88:93,95:96,98:99,104,108:114,122:124,126:127,129)] <- "SDunes"
#allpops[c(102:103,105:107,131:134,136:138)] <- "NCali"
#allpops[c(115:116,135)] <- "SCoast"
#allpops[c(121)] <- "SOreg"
cbind(indNames(snpgenind), allpops)

## calculate AFs for each pop...let's just run for pops with N>=5
#ncoast_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="NCoast")+1)])[,4]
ndunes_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="N Dunes")+1)])[,4]
sdunes_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="S Dunes")+1)])[,4]
#scoast_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="S Coast")+1)])[,4]
soreg_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="S Oregon")+1)])[,4]
ncali_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="N California")+1)])[,4]
orcasc_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="ORcasc")+1)])[,4]
bmt_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="ORbmt")+1)])[,4]
lssn_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="Lssn")+1)])[,4]
co_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)=="CO")+1)])[,4]

coma_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)%in%c("NCoast","NDunes","SDunes","SCoast","SOreg","NCali","Trinidad"))+1)])[,4]
mtn_AF <- maf(snpvcf[,c(1,which(pop(snpgenind)%in%c("ORcasc","WAcasc","WY","ORbmt","Lssn","CO")))])
all_AF <- maf(snpvcf)[,4]


## Step 1: simulate source population. Let's just run one source population at a time...


## Step 2: simulate alleles for individual based on frequencies observed in population
# let's use allele freq info from all genotyped individuals in the GBAS panel for better results
sim_genos <- function(AF, N=1000, sumAF=all_AF){
  simgeno <- as.data.frame(matrix(NA,nrow=N,ncol=length(AF)))
  names(simgeno) <- names(AF)
  for (snp in 1:length(AF)){
    maf <- AF[snp]
    if(is.na(maf)) maf <- sumAF[snp]
    rand1 <- runif(n=N, 0,1)
    rand2 <- runif(n=N, 0,1)
    simgeno[which(rand1<maf & rand2<maf), snp] <- 0
    simgeno[which(rand1<maf & rand2>=maf), snp] <- 0.5
    simgeno[which(rand1>=maf & rand2<maf), snp] <- 0.5
    simgeno[which(rand1>=maf & rand2>=maf), snp] <- 1
  }#for
  return(simgeno)
}#sim_genos function


## Step 3: calculate fORCA(Sm)
calc_fORCA <- function(simgeno, AF, K){
  forca <- c() #set up empty list
  
  for (i in 1:nrow(simgeno)){
    p <- c()# set up empty list...
    # calculate genotype probability for each marker
    for (m in 1:ncol(simgeno)){
      # calculate probability of genotype based on marker AFs in each pop
      geno <- as.numeric(simgeno[i,m])
      maf <- AF[m]
      if(is.na(maf)) maf <- all_AF[m]
      if(maf==0) maf<-0.01
      if(geno==0){
        pm <- (1/K)*2*maf
      }#if
      if (geno==0.5){
        pm <- (1/K)*(1)*maf*(1-maf)
      }#if
      if (geno==1){
        pm <- (1/K)*(2*(1-maf))
      }#if
      p <- append(p, pm)
    }#for markers
    
    # calculate overall genotype probability
    forca_i <- prod(p)
    forca <- append(forca, forca_i)
  }#for individuals
  return(forca)
}#calc_fORCA function


run_fORCA_allpops <- function(simgeno, pop){
  probs <- data.frame(id=rownames(simgeno), pop=pop, NDunes=NA, SDunes=NA, SOreg=NA, NCali=NA, Border=NA, 
                      ORcasc=NA, Lassen=NA, BMT=NA, CO=NA, Humboldt=NA, Mountain=NA)
  probs$NDunes <- calc_fORCA(simgeno, ndunes_AF, K=9)
  probs$SDunes <- calc_fORCA(simgeno, sdunes_AF, K=9)
  probs$SOreg <- calc_fORCA(simgeno, soreg_AF, K=9)
  probs$NCali <- calc_fORCA(simgeno, ncali_AF, K=9)
  probs$Border <- calc_fORCA(simgeno, border_AF, K=9)
  probs$ORcasc <- calc_fORCA(simgeno, orcasc_AF, K=9)
  probs$Lassen <- calc_fORCA(simgeno, lssn_AF, K=9)
  probs$BMT <- calc_fORCA(simgeno, bmt_AF, K=9)
  probs$CO <- calc_fORCA(simgeno, co_AF, K=9)
  probs$Humboldt <- calc_fORCA(simgeno, coma_AF, K=2)
  probs$Mountain <- calc_fORCA(simgeno, mtn_AF, K=2)
  return(probs)
}#run_fORCA_allpops


## Simulate each population and calculate fORCAs
probs <- data.frame(id=rep(1:1000, 9), pop=NA, NDunes=NA, SDunes=NA, SOreg=NA, NCali=NA, Border=NA,
                    ORcasc=NA, Lassen=NA, BMT=NA, CO=NA, Humboldt=NA, Mountain=NA)

ndunes_sim <- sim_genos(ndunes_maf$MAF)
sdunes_sim <- sim_genos(sdunes_maf$MAF)
soreg_sim <- sim_genos(soreg_maf$MAF)
ncali_sim <- sim_genos(ncali_maf$MAF)
border_sim <- sim_genos(border_maf$MAF)
orcasc_sim <- sim_genos(orcasc_maf$MAF)
bmt_sim <- sim_genos(bmt_maf$MAF)
lssn_sim <- sim_genos(lssn_maf$MAF)
co_sim <- sim_genos(co_maf$MAF)

probs[1:1000,] <- run_fORCA_allpops(ndunes_sim, "NDunes")
probs[1001:2000,] <- run_fORCA_allpops(sdunes_sim, "SDunes")
probs[2001:3000,] <- run_fORCA_allpops(soreg_sim, "SOreg")
probs[3001:4000,] <- run_fORCA_allpops(ncali_sim, "NCali")
probs[4001:5000,] <- run_fORCA_allpops(border_sim, "Border")
probs[5001:6000,] <- run_fORCA_allpops(orcasc_sim, "ORCasc")
probs[6001:7000,] <- run_fORCA_allpops(lssn_sim, "Lassen")
probs[7001:8000,] <- run_fORCA_allpops(bmt_sim, "BMT")
probs[8001:9000,] <- run_fORCA_allpops(co_sim, "CO")


## assign each individual to a pop and subspecies based on highest probability
probs$Assigned <- NA
probs$Subspecies <- NA
for (i in 1:nrow(probs)){
  pops <- names(probs)[3:10]
  maxpop <- which(probs[i,3:11]==max(probs[i,3:11]))
  probs$Assigned[i] <- sample(pops[maxpop], size=1)
  if(probs$Humboldt[i] > probs$Mountain[i]) probs$Subspecies[i] <- "Humboldt"
  if(probs$Humboldt[i] < probs$Mountain[i]) probs$Subspecies[i] <- "Mountain"
  if(probs$Humboldt[i]==probs$Mountain[i]) probs$Subspecies[i] <- sample(c("Humboldt","Mountain"), size=1)
}#for

##
table(probs$pop, probs$Assigned)
table(probs$pop, probs$Subspecies)



