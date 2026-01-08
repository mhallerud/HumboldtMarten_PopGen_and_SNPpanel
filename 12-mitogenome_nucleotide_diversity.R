####################################################
### MITOGENOME NUCLEOTIDE DIVERSITY CALCULATIONS ###
####################################################
# load dependencies
library(pegas)
library(adegenet)
library(ggplot2)

# read in manually curated table with sampleID, population info/citation/etc 
# associated with GenBank accession for each mitogenome seq.
info <- read.csv("maca_mitogenome_info.csv")

############################################
### FULL MITOGENOME NUCLEOTIDE DIVERSITY ###
############################################
## read in raw mitogenomes
# NOTE: These were first aligned and trimmed via MAFFT algorithm in Geneious
mitoseqs <- adegenet::fasta2DNAbin('AlignedTrimmed_MartenMitogenomes.fasta')
ids <- unlist(lapply(rownames(mitoseqs), function(x) strsplit(x," ")[[1]][1]))#extract IDs

# grab matching pop info
pops <- info$Population[match(ids,info$SequenceID)]

# remove partial mitogenome sequences
type <- info$Type[match(ids,info$SequenceID)]
mitoseqs <- mitoseqs[-which(type=="partial"),]
popsub <- pops[-which(type=="partial")]
mitoseqs#102 remaining


## subset sequences per pop/region (only if N>=2)
table(popsub)
qci <- mitoseqs[which(popsub=="BC- Queen Charlotte Islands"),]
bmt <- mitoseqs[which(popsub=="Blue Mountains"),]
lassen <- mitoseqs[which(popsub=="CA - Lassen"),]
shasta <- mitoseqs[which(popsub=="Califorina - Shasta"),]
cacoast <- mitoseqs[which(popsub=="California Coast"),]
sierra <- mitoseqs[which(popsub=="California Sierra"),]
orcasc <- mitoseqs[which(popsub=="Oregon Cascades"),] 
orcoast <- mitoseqs[which(popsub=="Oregon Coast"),]
wacasc <- mitoseqs[which(popsub=="Washington- Cascades"),]
olymp <- mitoseqs[which(popsub=="Washington- Olympia"),]

# calculate mitogenome-level pi for each pop
pi_orcoast <- pegas::nuc.div(orcoast, variance=TRUE)#7.884e-05
pi_cacoast <- pegas::nuc.div(cacoast, variance=TRUE)#1.410e-04
pi_orcasc <- pegas::nuc.div(orcasc, variance=TRUE)#3.707e-04
pi_shasta <- pegas::nuc.div(shasta, variance=TRUE)#5.353e-05
pi_sierra <- pegas::nuc.div(sierra, variance=TRUE)#7.45e-04
pi_wacasc <- pegas::nuc.div(wacasc, variance=TRUE)#1.82e-03
pi_olympia <- pegas::nuc.div(olymp, variance=TRUE)#1.51e-04
pi_lassen <- pegas::nuc.div(lassen, variance=TRUE)
pi_bmt <- pegas::nuc.div(bmt, variance=TRUE)
pi_qci <- pegas::nuc.div(qci, variance=TRUE)

# calculate overall pis
humboldt <- mitoseqs[which(popsub%in%c("Oregon Coast","California Coast")),]
wmontane <- mitoseqs[which(popsub%in%c("Oregon Cascades","CA - Lassen")),]
pi_humb <- pegas::nuc.div(humboldt, variance=TRUE)
pi_wmontane <- pegas::nuc.div(wmontane, variance=TRUE)
pi_all <- pegas::nuc.div(mitoseqs, variance=TRUE)

# combine into table
pi <- rbind(pi_orcoast, pi_cacoast, pi_orcasc, pi_shasta, pi_sierra, pi_wacasc, 
            pi_olympia, pi_bmt, pi_lassen, pi_qci, pi_humb, pi_wmontane, pi_all)
pi <- as.data.frame(pi)
names(pi) <- c("Pi_Est","Pi_Var")
pi$Population <- c("OR Coast","CA Coast","OR Cascades","CA Shasta","CA Sierra","WA Cascades",
                   "Olympia", "Blue Mountains", "Lassen","Haida Gwaii","Humboldt","W Montane","Pacific Martens")
pi$N <- c(nrow(orcoast), nrow(cacoast), nrow(orcasc), nrow(shasta), nrow(sierra), nrow(wacasc), 
                nrow(olymp), nrow(bmt), nrow(lassen), nrow(qci), nrow(humboldt), nrow(wmontane), nrow(maca_coi))



#######################################
### COI REGION NUCLEOTIDE DIVERSITY ###
#######################################
# read in data and grab ids + po
# these were extracted from all MAFFT-aligned sequences in Geneious
maca_coi <- adegenet::fasta2DNAbin('maca_mitogenomes_COIregion.fasta')
maca_coi#122 samples x 493 bp
coi_ids <- unlist(lapply(rownames(maca_coi), function(x) strsplit(x," ")[[1]][1]))
coi_pops <- info$Population[match(coi_ids,info$SequenceID)]

# subset per pop/region
table(coi_pops)
orcoast_coi <- maca_coi[which(coi_pops=="Oregon Coast"),]
cacoast_coi <- maca_coi[which(coi_pops=="California Coast"),]
orcasc_coi <- maca_coi[which(coi_pops=="Oregon Cascades"),] 
wacasc_coi <- maca_coi[which(coi_pops=="Washington- Cascades"),]
shasta_coi <- maca_coi[which(coi_pops=="Califorina - Shasta"),]
sierra_coi <- maca_coi[which(coi_pops=="California Sierra"),]
olymp_coi <- maca_coi[which(coi_pops=="Washington- Olympia"),]
bmt_coi <- maca_coi[which(coi_pops=="Blue Mountains"),]
lassen_coi <- maca_coi[which(coi_pops=="CA - Lassen"),]
montana_coi <- maca_coi[which(coi_pops=="-"),]
bc_coi <- maca_coi[which(coi_pops=="BC- Queen Charlotte Islands"),]

# calculate COI region pi for each region
coipi_orcoast <- pegas::nuc.div(orcoast_coi, variance=TRUE)
coipi_cacoast <- pegas::nuc.div(cacoast_coi, variance=TRUE)
coipi_orcasc <- pegas::nuc.div(orcasc_coi, variance=TRUE)
coipi_wacasc <- pegas::nuc.div(wacasc_coi, variance=TRUE)
coipi_shasta <- pegas::nuc.div(shasta_coi, variance=TRUE)
coipi_sierra <- pegas::nuc.div(sierra_coi, variance=TRUE)
coipi_olymp <- pegas::nuc.div(olymp_coi, variance=TRUE)
coipi_bmt <- pegas::nuc.div(bmt_coi, variance=TRUE)
coipi_lassen <- pegas::nuc.div(lassen_coi, variance=TRUE)
coipi_bc <- pegas::nuc.div(bc_coi, variance=TRUE)
coipi_montan <- pegas::nuc.div(montana_coi, variance=TRUE)

# calculate overall pis
humb_coi <- maca_coi[which(coi_pops%in%c("Oregon Coast","California Coast")),]
wmont_coi <- maca_coi[which(coi_pops%in%c("CA - Lassen","Oregon Cascades")),]
coipi_humb <- pegas::nuc.div(humb_coi, variance=TRUE)
coipi_wmont <- pegas::nuc.div(wmont_coi, variance=TRUE)
coipi_all <- pegas::nuc.div(maca_coi, variance=TRUE)

# convert to DF
coipi <- rbind(coipi_orcoast, coipi_cacoast, coipi_orcasc, coipi_wacasc, coipi_shasta, 
               coipi_sierra, coipi_olymp, coipi_bmt, coipi_lassen, coipi_bc, coipi_montan, coipi_humb, coipi_wmont, coipi_all)
coipi <- as.data.frame(coipi)
names(coipi) <- c("COI_Pi_Est", "COI_Variance")
coipi$Population <- c("OR Coast","CA Coast","OR Cascades","WA Cascades","CA Shasta","CA Sierra","Olympia",
                      "Blue Mountains","Lassen","Haida Gwaii","Montana hybrid zone", "Humboldt","W Montane","Pacific Martens")
coipi$COI_N <- c(nrow(orcoast_coi), nrow(cacoast_coi), nrow(orcasc_coi), nrow(wacasc_coi), 
                 nrow(shasta_coi), nrow(sierra_coi), nrow(olymp_coi), nrow(bmt_coi), nrow(lassen_coi), nrow(bc_coi), nrow(montana_coi),
                 nrow(humb_coi), nrow(wmont_coi), nrow(maca_coi))

# merge mitogenome and COI pi estimates into single DF
pi_merged <- merge(pi, coipi, by="Population", all=TRUE)
View(pi_merged)
write.csv(pi_merged, 'MACA_mitogenome_pi_estimates.csv')



##################
### FINAL PLOT ###
##################
pi <- mitopi
pi$Citation <- c(rep("Schwartz et al. 2020",9),"Colella et al. 2021",
                     rep("Schwartz et al. 2020",2),"Schwartz et al. 2020; Colella et al. 2021")
pi$Group <- c("Humboldt pop","Humboldt pop","Montane pop",rep(NA,7),"Humboldt",NA,"Pacific Marten")
pi$Species <- c("M.c. humboldtensis (OR)","M.c. humboldtensis (CA)", rep(NA,8),
                "M.c.humboldtensis",NA,"Martes caurina")
pi <- pi[,c("Population","Species","Group","Pi_Est","Pi_Var","N","Citation")]

# add values from other species mined from the literature
allpi <- rbind(pi,
      c("Striped hyena","Hyaena hyaena","VU",2.3e-4,NA,14,"Westbury et al. 2018"),#https://academic.oup.com/mbe/article/35/5/1225/4924857#118216583
      c("Island gray fox (Santa Cruz)","Urocyon littoralis (S. Cruz)","NT",2.4e-4, NA, 42, "Hofman et al. 2015"),#https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0118240&type=printable (all island gray foxes=0.00113)
      c("European mink (Fr-Sp)","Mustela lutreola","CR",2.4e-4, (3e-5)^2, 42, "Skorupski et al. 2024"),#https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/csp2.13291
      c("Black-footed ferret","Mustela nigripes","EN",2.99e-4,NA,4,"Etherington et al. 2022"),#https://academic.oup.com/jhered/article/113/5/500/6657699
      c("Iberian lynx","Lynx pardinus","VU",3.4e-4, (3.2e-4)^2, 59,"Casas-Marce et al. 2013"),#https://doi-org.oregonstate.idm.oclc.org/10.1111%2Fmec.12498
      c("Malayan Tiger","Panthera tigris jacksoni","EN",32/193 * 2.47e-3,NA,9,"Liu et al. 2018"),#https://www.cell.com/current-biology/fulltext/S0960-98221831214-4#mmc3
      c("Mediterranean monk seal","Monachus monachus","EN",5e-4,NA,14,"Rey-Iglesia et al. 2021"),#https://academic.oup.com/zoolinnean/article/191/4/1147/5897344
      c("Cheetah (Namibia)","Acinonyx jubatus (Namibia)","VU", 7.1e-4,(3.3e-4)^2, 4,"Dobrynin et al. 2015"),#https://link.springer.com/article/10.1186/s13059-015-0837-4#MOESM2
      c("Iberian Wolf", "Canis lupus signatus", "NT",8.2e-4, (4.3e-4)^2, 29, "Salado et al. 2023"),#https://www.mdpi.com/2073-4425/14/1/75
      c("Fisher","Pekania pennanti","LC",8.8e-4, NA, 40, "Knaus et al. 2011"),#https://link.springer.com/article/10.1186/1472-6785-11-10
      c("African wild dogs","Lycaon pictus","EN",9e-4,NA,20,"Tensen et al. 2025"),#https://www.tandfonline.com/doi/full/10.1080/24701394.2025.2558612#abstract
      c("Himalyan red panda","Ailurus styani","EN",1.005e-3,NA,13,"Hu et al. 2020"),#https://www.science.org/doi/full/10.1126/sciadv.aax5751#F2
      c("Eurasian otter (UK)","Lutra lutra (UK)","NT",1.6e-3,NA,46,"Du Plessis et al. 2025"),#https://link.springer.com/article/10.1007/s12686-025-01383-9
      c("Snow leopard","Panthera uncia","VU",1.77e-3, (0.00008)^2, 136,"Wang et al. 2025"),#https://link-springer-com.oregonstate.idm.oclc.org/content/pdf/10.1007/s10592-024-01658-y.pdf
      c("Chinese red panda","Ailurus fulgens","EN",3.7e-3,NA,36,"Hu et al. 2020"),#https://www.science.org/doi/full/10.1126/sciadv.aax5751#F2
      c("Flat-headed cat","Prionailurus planiceps","EN",3.96e-3,(4.8e-4)^2,21,"Patel et al. 2017"),#https://link-springer-com.oregonstate.idm.oclc.org/article/10.1007/s10592-017-0990-2
      c("Siberian weasel","Mustela siberica","LC",5.19e-3,NA,20,"Shalabi 2016"),#https://eprints.lib.hokudai.ac.jp/repo/huscap/all/67137/MOHAMMED_AMIN_MOHAMMED_MOHAMMED_SHALABI.pdf
      c("Sable","Martes zibellina","LC",5.8e-3,(1e-4)^2,75,"Li et al. 2021"),#https://link-springer-com.oregonstate.idm.oclc.org/content/pdf/10.1007/s42991-020-00092-0.pdf
      c("Japanese weasel","Mustela itatsi","NT",6.19e-3,NA,26,"Shalabi 2016")#https://eprints.lib.hokudai.ac.jp/repo/huscap/all/67137/MOHAMMED_AMIN_MOHAMMED_MOHAMMED_SHALABI.pdf
)

# prep fields for plotting...
groups <- allpi$Group
groups[groups%in%c("VU","NT","LC","EN","CR")] <- "Other"
groups <- factor(groups, levels=c("Humboldt pop","Humboldt","Pacific Marten","Other"))

allpi$Pi_Est <- as.numeric(allpi$Pi_Est)
allpi$Pi_Var <- as.numeric(allpi$Pi_Var)

allpi$Population <- factor(allpi$Population,levels=allpi$Population)
allpi$Species <- factor(allpi$Species, levels=allpi$Species)

# plot
ggplot(allpi[-which(is.na(groups)),])+
  geom_bar(aes(x=Species, y=Pi_Est, fill=groups[-which(is.na(groups))]), stat='identity')+
  theme_classic()+
  scale_fill_manual(values=c('plum','purple2',"purple4","gray"))+
  ggtitle(expression(paste("Mitogenome Nucleotide Diversity (", pi,")")))+
  ylab("mean pairwise differences")+
  scale_y_continuous(expand=c(0,0))+
  xlab(NULL)+theme(legend.position="none", text=element_text(size=14),
                   axis.text.y=element_text(face="italic"))+
                   #axis.text.x=element_text(angle=90,hjust=1))+
  coord_flip()
#save as 600x400



