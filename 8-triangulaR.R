##triangulaR documentation: https://github.com/omys-omics/triangulaR
library(triangulaR)
library(vcfR)
library(stringr)

# STEP 1: read in data
#sampleinfo <- read.csv("../SampleLocations/Marten_Hets_Fhom_coords.csv", stringsAsFactors=FALSE, fileEncoding="Latin1")
vcf <- read.vcfR("LDprunedRAD.vcf")
plate <- c(rep("Plate1",27), rep("Plate2",16), rep("Plate3",21))


# STEP 2: set up popmap
# grab sample IDs and popnames from VCF
sampleids <- colnames(vcf@gt)[-1]
popmap <- data.frame(id=sampleids, pop=str_split(sampleids,"_",n=2,simplify=TRUE)[,1])
unique(popmap$pop)#check

# rename pops (only those that will be labeled)
popmap$pop[popmap$pop=="NCali"] <- "N California"
popmap$pop[popmap$pop=="NDunes"] <- "N Dunes"
popmap$pop[popmap$pop=="ORcasc"] <- "OR Cascades"
popmap$pop[popmap$pop=="SDunes"] <- "S Dunes"
popmap$pop[popmap$pop=="SOreg"] <- "S Oregon"
popmap$pop[popmap$pop=="Trinidad"] <- "N California"
popmap$pop[popmap$pop=="Lssn"] <- "Lassen"

# set up e. montane vs. humboldt popmap
subspecies <- popmap
subspecies$pop[subspecies$pop%in%c("S Oregon","N California","Border","MapleCrk","S Dunes","N Dunes")] <- "Humboldt"
subspecies$pop[subspecies$pop%in%c("WAcasc","ORbmt","WY","CO")] <- "E. Montane"
unique(subspecies$pop)#pop



# STEP 3: choose sites above a given allele frequency difference threshold
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
# it's recommended to test different thresholds here

# Humboldt vs eastern Montane- 183 sites
diff90_hm <- alleleFreqDiff(vcfR=vcf, pm=subspecies, p1="Humboldt", p2="E. Montane", difference=0.9)#7
diff80_hm <- alleleFreqDiff(vcfR=vcf, pm=subspecies, p1="Humboldt", p2="E. Montane", difference=0.8)#52
diff70_hm <- alleleFreqDiff(vcfR=vcf, pm=subspecies, p1="Humboldt", p2="E. Montane", difference=0.7)#180


# run for Humboldt vs. OR Cascades
humbcasc <- vcf[,c(1,which(subspecies$pop%in%c("OR Cascades","Humboldt"))+1)]
hcpop <- subspecies[subspecies$pop%in%c("OR Cascades","Humboldt"),]
diff90_hc <- alleleFreqDiff(vcfR=humbcasc, pm=hcpop, p1="Humboldt", p2="OR Cascades", difference=0.9)#22
diff80_hc <- alleleFreqDiff(vcfR=humbcasc, pm=hcpop, p1="Humboldt", p2="OR Cascades", difference=0.8)#66
diff70_hc <- alleleFreqDiff(vcfR=humbcasc, pm=hcpop, p1="Humboldt", p2="OR Cascades", difference=0.7)#197



# STEP 4: Calculate hybrid index and heterozygosity for each sample. Values are returned in a data.frame
hi.het90_hm <- hybridIndex(vcfR=diff90_hm, pm=subspecies, p1="Humboldt", p2="E. Montane")
hi.het80_hm <- hybridIndex(vcfR=diff80_hm, pm=subspecies, p1="Humboldt", p2="E. Montane")
hi.het70_hm <- hybridIndex(vcfR=diff70_hm, pm=subspecies, p1="Humboldt", p2="E. Montane")

hi.het90_casc <- hybridIndex(vcfR=diff90_hc, pm=hcpop, p1="Humboldt", p2="OR Cascades")
hi.het80_casc <- hybridIndex(vcfR=diff80_hc, pm=hcpop, p1="Humboldt", p2="OR Cascades")
hi.het70_casc <- hybridIndex(vcfR=diff70_hc, pm=hcpop, p1="Humboldt", p2="OR Cascades")



# STEP 5: Visualize as triangle plot
# compare deltas for Humboldt vs. eastern montane
triangle.plot(hi.het90_hm)#delta=0.90
triangle.plot(hi.het80_hm)#delta=0.80
triangle.plot(hi.het70_hm)#delta=0.70
#check effect of missingness
missing.plot(hi.het70_hm)
# final plot using delta 0.70
hi.het70_hm$pop <- factor(hi.het70_hm$pop, levels=c("Humboldt","E. Montane","OR Cascades","Lassen"))
triangle.plot(hi.het70_hm, alpha=0.7, colors=c("blue","red","turquoise","purple4"), cex=3)+
  theme(legend.position="left", text=element_text(size=14), legend.title=element_blank())
#save 500x300

# check delta for Humboldt vs. OR Cascades
triangle.plot(hi.het90_casc)
triangle.plot(hi.het80_casc)
triangle.plot(hi.het70_casc)
#check effect of missingness
missing.plot(hi.het70_casc)
# final plot using delta 0.70
triangle.plot(hi.het70_casc, alpha=0.7, colors=c("blue","turquoise"), cex=3, ind.labels=TRUE)+
  theme(legend.position="right", text=element_text(size=14), legend.title=element_blank())
#save 500x300


### INTERPRETATION:
#https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.14039
#IBD: straight line between two lower clusters
#Neutral diffusion (i.e. gene flow): samples follow curve
#No introgression: only F1 hybrids at ~0.5/0.5 (inflection point on curve) 


