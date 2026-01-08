#### GBAS SNP PANEL ASSESSMENT ####
### PURPOSE: Generate GT-seq plots for allele balances, per-locus and per-sample OTRs, 
####      genotyping rates, and PID metrics
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
options(warn=0) #avoid ggplot warnings being converted into errors



##########################
### REFORMAT GENOTYPES ###
##########################
# read in data
genos <- read.csv("GBAS_genotypes.csv")

# extract consensus genotypes
consensus <- genos[which(genos$Filepath=="CONSENSUS"),]
dim(consensus)#480 consensus genotypes
unique(consensus$Type)
names(consensus)

# clean up genotypes- only MACA polymorphic sites (no SRY)
genos <- consensus[,c(14:106,109:117)]
unique(unlist(genos))
genos[genos==0] <- NA
genos[genos==""] <- NA
genos[genos=="FALSE"] <- NA
genos[genos=="9"] <- NA
genos[genos=="O"] <- NA
genos[genos=="NA"] <- NA
genos[genos==""] <- NA
genos[genos=="-"] <- NA
# these are for SRY....
genos[genos=="GG* (1 rep)"] <- "GG"
genos[genos=="GG (1 rep)"] <- "GG"
genos[genos=="GG (2 reps)"] <- "GG"
genos[genos=="GG (2/6 reps)"] <- "GG"
genos[genos=="AA* (1 rep)"] <- "AA"
genos[genos=="GG "] <- "GG"
sort(unique(unlist(genos)))
consensus <- cbind(consensus[,1:14], genos)
# remove duplicate SNP field - this SNP was genotyped twice with more specific probe
# given amplification of Myodes
consensus <- consensus[,-which(names(consensus)=="MACA_67356.0.111.1")]

# calculate genotyping rate (per-amplicon so only "main" SNPs)
consensus$N_GT <- apply(consensus[,15:108], 1, FUN=function(x) sum(!is.na(x)))
hist(consensus$N_GT)

# classify into successful, marginal/rerun, and failed scats based on # SNPs and # reps
consensus$Failed[consensus$N_GT>=75] <- "Successful" #80% genotyped
#consensus$Failed[which(consensus$N_GT<75 & consensus$N_GT>=47)] <- "Marginal" #50-80% SNPs genotyped
consensus$Failed[consensus$N_GT<75] <- "Failed" #<80% genotyped
consensus$Failed[which(consensus$N_GT>=75 & consensus$Rep>3)] <- "Rerun"
consensus$Failed <- factor(consensus$Failed, levels=c("Failed","Rerun","Marginal","Successful"))

# clean up sample types
unique(consensus$Type)
consensus$Type[consensus$Type=="hair_marten"] <- "known_scat" #only 2 hair samples, bin with scats

# remove duplicate samples (genotyped on 2 plates)
consensus <- consensus[consensus$PlateID!="MACA_P10B",] # P10 included twice in CSV -remove
tbl <- as.data.frame(table(unique(consensus$SampleID)))
tbl[which(tbl$Freq>1),]
consensus$SampleID[consensus$SampleID=="KMSM83"] <-  c("KMSM83A","KMSM83B")
consensus$SampleID[consensus$SampleID=="KMSM82"] <- c("KMSM82A","KMSM82B")

# summary of genotyping types
dim(consensus)#431
table(consensus$Type)#90 tissue, 258 scat, 53 known scat, 18 blank, 12 nontarget species

# summary of genotyping success
# all marten samples - 71% success + 3.5% rerun
table(consensus$Failed[consensus$Type%in%c("scat_marten","tissue_marten","known_scat")]) / sum(consensus$Type%in%c("scat_marten","tissue_marten","known_scat"))
# tissue only - 90% success + 0% rerun
table(consensus$Failed[consensus$Type=="tissue_marten"]) / sum(consensus$Type=="tissue_marten")
# all scats - 66% success + 4.5% rerun
table(consensus$Failed[consensus$Type%in%c("scat_marten","known_scat")]) / sum(consensus$Type%in%c("scat_marten","known_scat"))
 # fresh scats only - 75% success + 4% rerun success
table(consensus$Failed[consensus$Type=="known_scat"]) / sum(consensus$Type=="known_scat")
# field-collected scats only - 64% success + 5% rerun success
table(consensus$Failed[consensus$Type=="scat_marten"]) / sum(consensus$Type=="scat_marten")

# clean up other fields...
consensus$On.Target.Reads <- as.numeric(consensus$On.Target.Reads)
consensus$Raw.Reads <- as.numeric(consensus$Raw.Reads)
consensus$X.GT <- consensus$N_GT / 92
consensus$X.On.Target <- consensus$On.Target.Reads / consensus$Raw.Reads
names(consensus)[c(3,6,8:11)] <- c("SeqRun","Reps","RawReads","OnTargetReads","PercentOnTarget","PercentGT")

# remove blanks, non-marten samples
consensus <- consensus[!consensus$Type%in%c("nontarget","blank"),]

# check
dim(consensus)#401
View(consensus)



#######################################
####### PLOTS FOR GTSEQ EVALUATION ####
#######################################
consensus$Type <- factor(consensus$Type, levels=c("scat_marten","known_scat","tissue_marten"))
blues <- rev(RColorBrewer::brewer.pal(3,"Blues"))

## % ON-TARGET READS BY SAMPLE TYPE AND FAILURE RATE
consensus$Type <- factor(consensus$Type, levels=c("tissue_marten","known_scat","scat_marten"))
p5 <- ggplot(data=consensus)+
  geom_boxplot(aes(y=PercentOnTarget, x=Failed, fill=Type))+
  theme_bw()+
  xlab(element_blank())+
  ylab(element_blank())+#"% On-Target Reads")+
  scale_y_continuous(labels=scales::percent)+
  scale_fill_manual(values=rev(blues),
                    breaks=c("tissue_marten","known_scat","scat_marten"),
                    labels=c("Tissue","Fresh Scat","Field Scat"))+
  theme(panel.grid.minor=element_blank(), text=element_text(size=12))+
  ggtitle("C")
#ggtitle("On-Target Sequencing Rate")


## ON-TARGET READS BY SAMPLE TYPE AND FAILURE RATE
p6 <- ggplot(data=consensus)+
  geom_boxplot(aes(y=OnTargetReads, x=Failed, fill=Type))+
  theme_bw()+
  ylab("On-target reads")+
  xlab("Sample Quality")+
  scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  scale_fill_manual(values=rev(blues),
                    breaks=c("tissue_marten","known_scat","scat_marten"),
                    labels=c("Tissue","Fresh Scat","Field Scat"))+
  theme(panel.grid.minor=element_blank(), text=element_text(size=12))+
  ggtitle("F")
#ggtitle("Reads by Sample Type")


## ON-TARGET READS BY SAMPLE TYPE
ggplot(data=consensus)+
  geom_boxplot(aes(y=OnTargetReads, x=Failed, fill=Type))+
  theme_bw()+
  ylab("Raw Read Counts")+
  xlab(element_blank())+
  scale_y_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+
  scale_fill_manual(values=blues,
                    breaks=c("tissue_marten","known_scat","scat_marten"),
                    labels=c("Tissue","Fresh Scat","Field Scat"))+
  theme(panel.grid.minor=element_blank(), text=element_text(size=18))+
  ggtitle("Reads by Sample Type")


## ON TARGET % PER SAMPLE
consensus$FailType <- paste(consensus$Failed, consensus$Type)
blues <- RColorBrewer::brewer.pal(3,"Blues")
#greens <- rev(RColorBrewer::brewer.pal(3, "Greens"))
#purples <- rev(RColorBrewer::brewer.pal(3, "Purples"))
consensus$SampleID <- factor(consensus$SampleID, levels=consensus$SampleID[order(consensus$PercentOnTarget)])
p3 <- ggplot(consensus)+#[consensus$Type%in%c("scat_marten","known_scat"),])+
  geom_bar(aes(x=SampleID, y=PercentOnTarget, col=Type, fill=Type),stat='identity')+
  scale_color_manual(values=blues, labels=c("Scat","Fresh Scat","Tissue"), breaks=c("scat_marten","known_scat","tissue_marten"))+
  scale_fill_manual(values=blues, labels=c("Scat","Fresh Scat","Tissue"), breaks=c("scat_marten","known_scat","tissue_marten"))+
  # scale_color_manual(values=c(greens, blues, purples),
  #                    labels=c("Successful (tissue)","Successful (fresh scat)","Successful (scat)",
  #                             "Rerun (fresh scat)","Rerun (scat)",
  #                             "Failed (tissue)",  "Failed (fresh scat)", "Failed (scat)"),
  #                    breaks=c("Successful tissue_marten","Successful known_scat","Successful scat_marten",
  #                             "Rerun known_scat","Rerun scat_marten",
  #                             "Failed tissue_marten","Failed known_scat","Failed scat_marten"))+
  # scale_fill_manual(values=c(greens,blues,purples),#skyblue
  #                   labels=c("Successful (tissue)","Successful (fresh scat)","Successful (scat)",
  #                            "Rerun (fresh scat)","Rerun (scat)",
  #                            "Failed (tissue)",  "Failed (fresh scat)", "Failed (scat)"),
#                   breaks=c("Successful tissue_marten","Successful known_scat","Successful scat_marten",
#                            "Rerun known_scat","Rerun scat_marten",
#                            "Failed tissue_marten","Failed known_scat","Failed scat_marten"))+
theme_bw()+
  scale_y_continuous(labels=scales::label_percent(), expand=c(0,0),limits=c(0,1))+
  scale_x_discrete(expand=c(0,0))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position='none',
        text=element_text(size=12), legend.title=element_blank(), 
        panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())+
  xlab(element_blank())+ylab(element_blank())+
  ggtitle("B")
#xlab("Sample")+ylab("On-Target Reads")#+
#ggtitle("Sample On-Target Rate")


## SNPs GENOTYPED PER SAMPLE
consensus$SampleID <- factor(consensus$SampleID, levels=consensus$SampleID[order(consensus$N_GT)])
p4 <- ggplot(consensus)+
  geom_bar(aes(x=SampleID, y=N_GT, col=Type, fill=Type),stat='identity')+
  scale_color_manual(values=RColorBrewer::brewer.pal(3,"Blues"), labels=c("Scat","Fresh Scat","Tissue"), breaks=c("scat_marten","known_scat","tissue_marten"))+
  scale_fill_manual(values=RColorBrewer::brewer.pal(3,"Blues"), labels=c("Scat","Fresh Scat","Tissue"), breaks=c("scat_marten","known_scat","tissue_marten"))+
  geom_hline(aes(yintercept=75), lty='dashed')+
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position='none',
        text=element_text(size=12), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())+
  xlab("Sample")+ylab("Loci genotyped")+
  ggtitle("E")
#ggtitle("Sample Genotyping Rate")


#### PER-LOCUS CALCULATIONS ####
gtseq <- read.csv("AllGenos.csv",
                  header=FALSE, 
                  col.names=c("RunID","Locus","A1","A2","Ratio","Geno","Call","Adj","Adj2","Total","OnTarget","IFI","Test","Notes","FailSeq","FailCount"))
head(gtseq)

# extract raw allele counts
gtseq$A1_reads <- as.numeric(unlist(lapply(strsplit(gtseq$A1,"="), "[", 2)))
gtseq$A2_reads <- as.numeric(unlist(lapply(strsplit(gtseq$A2,"="), "[", 2)))

## (RE) calculate IFI
gtseq$IFI[which(gtseq$Call=="HET")] <- NA
gtseq$IFI[which(gtseq$Call=="A1HOM")] <- gtseq$A2_reads[which(gtseq$Call=="A1HOM")] / 
  (gtseq$A1_reads[which(gtseq$Call=="A1HOM")]+gtseq$A2_reads[which(gtseq$Call=="A1HOM")])
gtseq$IFI[which(gtseq$Call=="A2HOM")] <- gtseq$A1_reads[which(gtseq$Call=="A2HOM")] / 
  (gtseq$A1_reads[which(gtseq$Call=="A2HOM")]+gtseq$A2_reads[which(gtseq$Call=="A2HOM")])

# add blank rows with file info
gtseq$Filename <- NA
gtseq$S_RawReads <- NA
gtseq$S_OnTargetReads <- NA
gtseq$S_OnTarget <- NA
gtseq$S_IFI <- NA

# find rows with file info
fileinfo <- gtseq[which(grepl("illumina",gtseq$Locus)),]

# add file info info into CSV
for (f in 1:nrow(fileinfo)){
  # ID which rows are associated with each file
  if(f!=nrow(fileinfo)){
    rows <- as.numeric(rownames(fileinfo)[f:(f+1)])
    rows[2] <- rows[2]-1
  }else{
    rows <- c(as.numeric(rownames(fileinfo)[f]), nrow(gtseq))
  }#ifelse
  
  # extract file info from header line
  filename <- fileinfo$Locus[f]
  rawreads <- fileinfo$A1[f]
  ontargetreads <- fileinfo$A2[f]
  percentontarget <- fileinfo$Ratio[f]
  ifi <- fileinfo$Geno[f]
  
  # add fileinfo to rows associated with sample 
  gtseq$Filename[rows[1]:rows[2]] <- filename
  gtseq$S_RawReads[rows[1]:rows[2]] <- strsplit(rawreads, ":")[[1]][2]
  gtseq$S_OnTargetReads[rows[1]:rows[2]] <- strsplit(ontargetreads,":")[[1]][2]
  gtseq$S_OnTarget[rows[1]:rows[2]] <- strsplit(percentontarget,":")[[1]][2]
  if (ifi==""){
    gtseq$S_IFI[rows[1]:rows[2]] <- mean(gtseq$IFI[rows[1]:rows[2]], na.rm=TRUE)
  }else{
    gtseq$S_IFI[rows[1]:rows[2]] <- strsplit(ifi,":")[[1]][2]
  }#ifelse
}#for
head(gtseq)
tail(gtseq)

## Correct sample name & PCR rep for these...
files <- read.csv("FilenameSampleID.csv")
gtseq$SampleID <- NA
gtseq$SampleID <- files$SampleID[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$PCRrep <- files$Rep[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$Plate <- files$PlateID[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$Run <- files$Run[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]

#clean up 
gtseq$S_RawReads <- as.numeric(gtseq$S_RawReads)
gtseq$S_OnTarget <- as.numeric(gtseq$S_OnTarget)
gtseq$S_OnTargetReads <- as.numeric(gtseq$S_OnTargetReads)
gtseq$S_IFI <- as.numeric(gtseq$S_IFI)

# convert filename info into sample #, indexes, sample ID, PCR rep info
basename <- basename(gtseq$Filename)
gtseq$SampleNumber <- unlist(lapply(strsplit(basename, "-"), "[", 2))
gtseq$i5 <- unlist(lapply(strsplit(basename,"-"), "[", 5))
gtseq$i7 <- unlist(lapply(strsplit(basename,"-"), "[", 6))
#gtseq$SampleID <- unlist(lapply(strsplit(basename,"_"), "[", 2))
gtseq$PCRrep <- unlist(lapply(strsplit(basename,"_"), "[", 3))

# remove rows with header info
fileinfo <- gtseq[which(grepl("fastq",gtseq$Locus)), c("Filename","SampleNumber","SampleID","PCRrep","S_RawReads","S_OnTargetReads","S_OnTarget","S_IFI")]
gtseq <- gtseq[-which(grepl("fastq",gtseq$Locus)),]

# clean up fields
gtseq$Geno <- factor(gtseq$Geno, levels=unique(gtseq$Geno))
gtseq$Call[which(is.na(gtseq$Call))] <- "00"

# remove duplicate loci (names vary per run...)
locusinfo <- read.csv("LocusInfo.csv")
locusinfo$FieldID[locusinfo$FieldID==""] <- NA
gtseqsub <- gtseq[which(gtseq$RunID=="Nextseq496" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq496"]),]
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq498" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq498"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq499" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq499"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq509" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq509"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq510" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq510"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq511" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq511"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq520" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq520"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq533" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq533"]),])
gtseqsub <- rbind(gtseqsub,
                  gtseq[which(gtseq$RunID=="Nextseq552" & gtseq$Locus%in%locusinfo$FieldID[locusinfo$Run=="Nextseq552"]),])
dim(gtseq); dim(gtseqsub)#check

# merge loci with different names
sort(unique(gtseqsub$Locus))
gtseqsub$Locus[gtseqsub$Locus=="SRY_MAAM.1"] <- "SRY_MAAM_1"
gtseqsub$Locus[gtseqsub$Locus=="SRY_MAAM.4"] <- "SRY_MAAM_4"
gtseqsub$Locus <- unlist(lapply(gtseqsub$Locus, function(x) strsplit(x,"\\.")[[1]][1]))
gtseqsub$Locus[gtseqsub$Locus=="MACA_67356:111"] <- "MACA_67356"
sort(unique(gtseqsub$Locus))

# check for loci with no reads
totreads_perlocus <- aggregate(data=gtseqsub, Total~Locus, FUN=sum)
noreads_loci <- totreads_perlocus$Locus[totreads_perlocus$Total==0]; noreads_loci
#gtseq <- gtseq[-which(gtseq$Locus %in% noreads_loci),]

# remove NTCs
gtseqsub <- gtseqsub[-which(startsWith(gtseqsub$SampleID,"NTC")),]

# calculate per-locus on-target rates
gtseq_agg <- aggregate(Total~Locus, data=gtseqsub[which(gtseqsub$SampleID%in%consensus$SampleID),], FUN=mean)
gtseq_agg$Locus <- factor(gtseq_agg$Locus, levels=gtseq_agg$Locus[order(gtseq_agg$Total)])

## read depth per locus
ggplot(gtseq)+
  geom_boxplot(aes(x=Locus, y=Total))+#, stat='identity')+
  theme_classic()+
  theme(axis.text.x=element_blank(), text=element_text(size=18))+
  ylab("Average read depth")+
  xlab("Locus")+
  ggtitle("Read Depths Across Loci")+
  ylim(0,5e4)

## average read depth per locus
ggplot(gtseq_agg[gtseq_agg$Total>10,])+
  geom_bar(aes(x=Locus, y=Total), stat='identity')+
  theme_classic()+
  theme(axis.text.x=element_blank(), text=element_text(size=18))+
  ylab("Average read depth")+
  xlab("Locus")+
  ggtitle("Read Depths Across Loci")

# average read depth per locus
ggplot(gtseqsub[gtseqsub$SampleID%in%consensus$SampleID,])+
  geom_boxplot(aes(x=Locus, y=Total/S_RawReads), col='red', alpha=0.5)+
  geom_boxplot(aes(x=Locus, y=Total*(OnTarget/100)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("Average Reads Per Locus")
#ylim(0, 2e05)


## PER-LOCUS ON-TARGET RATE (AVG)
# calculate average on-target rate per locus
gtseq_agg <- aggregate(OnTarget/100 ~ Locus, data=gtseqsub, FUN=mean)
names(gtseq_agg) <- c("Locus","mean")
gtseq_agg$Locus <- factor(gtseq_agg$Locus, levels=gtseq_agg$Locus[order(gtseq_agg$mean)])
gtseq_agg$Type <- NA
gtseq_agg$Type[gtseq_agg$Locus%in%c("SRY_MAAM_1","SRY_MAAM_4")] <- "SRY"
p1 <- ggplot(gtseq_agg)+
  geom_bar(aes(x=Locus, y=mean, fill=Type, col=Type), stat="identity")+
  theme_bw()+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), 
        axis.text.x=element_blank(), text=element_text(size=12), legend.position="none",
        axis.ticks.x=element_blank())+
  ylim(0,1)+
  scale_y_continuous(labels=scales::percent, expand=c(0,0))+
  ylab("On-Target Reads")+
  xlab(element_blank())+
  ggtitle("A")
#+ggtitle("Locus On-Target Rate")

ggplot(gtseqsub)+
  #geom_boxplot(aes(x=Locus, y=Total/S_RawReads), col='red', alpha=0.5)+
  geom_boxplot(aes(x=Locus, y=(OnTarget)/100))+
  geom_point(data=gtseq_agg, aes(x=Locus, y=mean), col='red', pch=8)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  ggtitle("Average % Reads per Locus")


### LOCUS GENOTYPING RATE
gtseq_agg <- aggregate(Call~Locus, gtseqsub, FUN=function(x) 1-(sum(x=="00")/length(x)))
gtseq_agg$Locus <- factor(gtseq_agg$Locus, levels=gtseq_agg$Locus[order(gtseq_agg$Call)])
gtseq_agg$Type <- NA
gtseq_agg$Type[gtseq_agg$Locus%in%c("SRY_MAAM_1","SRY_MAAM_4")] <- "SRY"
p2 <- ggplot(gtseq_agg)+
  geom_bar(aes(x=Locus,y=Call,col=Type, fill=Type), stat="identity")+
  theme_bw()+
  xlab("Locus")+
  scale_y_continuous(labels=scales::percent, expand=c(0,0), limits=c(0,.8))+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none", axis.text.x=element_blank(), text=element_text(size=12),
        axis.ticks.x=element_blank())+
  ylab("Samples genotyped")+
  ggtitle("D")

# on-target read breakdown by loci
unique(gtseq$Locus)# check for "failed" sequences showing up in loci...
ggplot(gtseqsub)+
  geom_boxplot(aes(x=Locus, y=Total/S_OnTargetReads), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  geom_hline(aes(yintercept=1/93), lty=2, col='red')+
  ggtitle("% Overall Reads per Locus")

# summary plot: average % ontarget per locus (order by ontarget)
gtseq_agg <- aggregate(OnTarget ~ Locus, data=gtseq, FUN=mean)
names(gtseq_agg) <- c("Locus","mean")
gtseq$Locus <- factor(gtseq$Locus, levels=gtseq_agg$Locus[order(gtseq_agg$mean)])

ggplot(gtseq[gtseq$SampleID%in%consensussub$SampleID,])+
  geom_boxplot(aes(x=Locus, y=OnTarget), fill='gray80')+
  #geom_point(data=gtseq_agg, aes(x=Locus, y=mean), col='blue', pch=8)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  #geom_hline(aes(yintercept=70), col='red', lty=2)+
  ggtitle("% OnTarget Per Locus")


## SET UP FULL FIGURE
fig1 <- ggarrange(p1,p3,p5,nrow=1,ncol=3,widths=c(1,1,1.2))
fig2 <- ggarrange(p2,p4,p6,nrow=1,ncol=3,widths=c(1,1,1.2))
#fig1 <- annotate_figure(fig1, top=text_grob("On-Target Rates", face = "bold", size = 14))
#fig2 <- annotate_figure(fig2, top=text_grob("Genotyping Rates", face="bold", size=14))
ggarrange(fig1, fig2, nrow=2, ncol=1)



##############################
#### ALLELE BALANCE PLOTS ####
##############################

#### INPUT: Concatenated *genos files ####
# read in "cat *genos" file
gtseq <- read.csv(#"/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2B_95pairs/July2024_Test/AllGenos_long.csv", 
                  "AllGenos.csv",
                  header=FALSE, 
                  col.names=c("Locus","A1","A2","Ratio","Geno","Call","Adj","Adj2","Total","OnTarget","IFI","Test","Notes","FailSeq","FailCount"))

#### STEP 0: CLEANING UP INPUT FILE ####
# extract raw allele counts
gtseq$A1_reads <- as.numeric(unlist(lapply(strsplit(gtseq$A1,"="), "[", 2)))
gtseq$A2_reads <- as.numeric(unlist(lapply(strsplit(gtseq$A2,"="), "[", 2)))
  
## (RE) calculate IFI
gtseq$IFI[which(gtseq$Call=="HET")] <- NA
gtseq$IFI[which(gtseq$Call=="A1HOM")] <- gtseq$A2_reads[which(gtseq$Call=="A1HOM")] / 
  (gtseq$A1_reads[which(gtseq$Call=="A1HOM")]+gtseq$A2_reads[which(gtseq$Call=="A1HOM")])
gtseq$IFI[which(gtseq$Call=="A2HOM")] <- gtseq$A1_reads[which(gtseq$Call=="A2HOM")] / 
  (gtseq$A1_reads[which(gtseq$Call=="A2HOM")]+gtseq$A2_reads[which(gtseq$Call=="A2HOM")])

# add blank rows with file info
gtseq$Filename <- NA
gtseq$S_RawReads <- NA
gtseq$S_OnTargetReads <- NA
gtseq$S_OnTarget <- NA
gtseq$S_IFI <- NA

# find rows with file info
fileinfo <- gtseq[which(grepl("fastq",gtseq$Locus)),]

# add file info info into CSV
for (f in 1:nrow(fileinfo)){
  # ID which rows are associated with each file
  if(f!=nrow(fileinfo)){
    rows <- as.numeric(rownames(fileinfo)[f:(f+1)])
    rows[2] <- rows[2]-1
  }else{
    rows <- c(as.numeric(rownames(fileinfo)[f]), nrow(gtseq))
  }#ifelse
  
  # extract file info from header line
  filename <- fileinfo$Locus[f]
  rawreads <- fileinfo$A1[f]
  ontargetreads <- fileinfo$A2[f]
  percentontarget <- fileinfo$Ratio[f]
  ifi <- fileinfo$Geno[f]
  
  # add fileinfo to rows associated with sample 
  gtseq$Filename[rows[1]:rows[2]] <- filename
  gtseq$S_RawReads[rows[1]:rows[2]] <- strsplit(rawreads, ":")[[1]][2]
  gtseq$S_OnTargetReads[rows[1]:rows[2]] <- strsplit(ontargetreads,":")[[1]][2]
  gtseq$S_OnTarget[rows[1]:rows[2]] <- strsplit(percentontarget,":")[[1]][2]
  if (ifi==""){
    gtseq$S_IFI[rows[1]:rows[2]] <- mean(gtseq$IFI[rows[1]:rows[2]], na.rm=TRUE)
  }else{
    gtseq$S_IFI[rows[1]:rows[2]] <- strsplit(ifi,":")[[1]][2]
  }#ifelse
}#for
head(gtseq)
tail(gtseq)

## Correct sample name & PCR rep for these...
files <- read.csv("FilenameSampleID.csv")
gtseq$SampleID <- NA
gtseq$SampleID <- files$SampleID[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$PCRrep <- files$Rep[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$Plate <- files$PlateID[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]
gtseq$Run <- files$Run[match(stringr::str_replace(basename(gtseq$Filename), ".fastq", ""), files$Filepath)]

#clean up 
gtseq$S_RawReads <- as.numeric(gtseq$S_RawReads)
gtseq$S_OnTarget <- as.numeric(gtseq$S_OnTarget)
gtseq$S_OnTargetReads <- as.numeric(gtseq$S_OnTargetReads)
gtseq$S_IFI <- as.numeric(gtseq$S_IFI)

# convert filename info into sample #, indexes, sample ID, PCR rep info
basename <- basename(gtseq$Filename)
gtseq$SampleNumber <- unlist(lapply(strsplit(basename, "-"), "[", 2))
gtseq$i5 <- unlist(lapply(strsplit(basename,"-"), "[", 5))
gtseq$i7 <- unlist(lapply(strsplit(basename,"-"), "[", 6))
#gtseq$SampleID <- unlist(lapply(strsplit(basename,"_"), "[", 2))
gtseq$PCRrep <- unlist(lapply(strsplit(basename,"_"), "[", 3))

# remove rows with header info
fileinfo <- gtseq[which(grepl("fastq",gtseq$Locus)), c("Filename","SampleNumber","SampleID","PCRrep","S_RawReads","S_OnTargetReads","S_OnTarget","S_IFI")]
gtseq <- gtseq[-which(grepl("fastq",gtseq$Locus)),]

# clean up fields
gtseq$Geno <- factor(gtseq$Geno, levels=unique(gtseq$Geno))
gtseq$Call[which(is.na(gtseq$Call))] <- "00"

# remove loci with no reads for clarity
totreads_perlocus <- aggregate(data=gtseq, Total~Locus, FUN=sum)
noreads_loci <- totreads_perlocus$Locus[totreads_perlocus$Total==0]; noreads_loci
gtseq <- gtseq[-which(gtseq$Locus %in% noreads_loci),]


#### STEP 1: Plot per-locus allele balance & read gtseq ####
plotAlleleBalance <- function(gtseq, locus, lim=NULL){
  sub <- gtseq[gtseq$Locus==locus,]
  if (sum(sub$Total)>0){
    #if((max(sub$A1_reads)>1000|max(sub$A2_reads)>1000)){
    #  xl<- yl<- 1000
    #}else{
    xl <- yl <- max(c(sub$A1_reads, sub$A2_reads))
    #}#ifelse
    ifelse(!is.null(lim), xl<-yl<- lim, xl<-yl<- max(c(sub$A1_reads, sub$A2_reads)))
    
    sub_het <- sub[which(sub$Call=="HET"),]
    sub_a1hom <- sub[which(sub$Call=="A1HOM"),]
    sub_a2hom <- sub[which(sub$Call=="A2HOM"&sub$A1_reads!=0),]
    het_slope <- mean(sub_het$A2_reads/sub_het$A1_reads)
    #a1hom_slope <- mean(sub_a1hom$A2_reads/sub_a1hom$A1_reads)
    #a2hom_slope <- mean(sub_a2hom$A2_reads/sub_a2hom$A1_reads)
    #a1hom_int <- quantile(sub_a1hom$A2_reads, 0.80, na.rm=T)
    #a2hom_int <- quantile(sub_a2hom$A1_reads, 0.80, na.rm=T)
    p <- ggplot()+
      geom_polygon(aes(x=c(0,xl,xl), y=c(0,0,yl/10)), fill='red', alpha=0.2)+
      geom_polygon(aes(x=c(0,0,xl/10), y=c(0,yl,yl)), fill='blue', alpha=0.2)+
      geom_polygon(aes(x=c(0,xl,xl,xl/5), y=c(0,yl/5,yl,yl)), fill='yellow', alpha=0.2)+
      geom_abline(aes(slope=1, intercept=0), lwd=0.1)+
      geom_point(data=sub, aes(x=A1_reads, y=A2_reads, col=Geno), alpha=0.3)+
      theme_classic()+
      ggtitle(locus)+
      xlab("A1")+ylab("A2")+
      theme(legend.position="none")+
      scale_x_continuous(expand=expansion(mult=0,add=0), limits=c(0,xl))+
      scale_y_continuous(expand=expansion(mult=0,add=0), limits=c(0,yl))
    #if (!is.na(a1hom_int)){
    #  p <- p + 
    #    geom_hline(aes(yintercept=a1hom_int), col="red", lwd=0.2, lty=2)
    #}#if
    #if (!is.na(a2hom_int)){
    #  p <- p+
    #  geom_vline(aes(xintercept=a2hom_int), col="blue",lwd=0.2, lty=2)
    #}#if
    #if (!is.na(het_slope)){
    #  p <- p+
    #    geom_abline(aes(slope=het_slope,intercept=0), lwd=1, col='yellow')
    #}#if
    #if (!is.na(a1hom_slope)){
    #  p <- p+
    #    geom_abline(aes(slope=a1hom_slope, intercept=0), lwd=0.2, col="red")
    #}#if
    #if (!is.na(a2hom_slope)){
    #  geom_abline(aes(slope=a2hom_slope, intercept=0.0), lwd=0.2, col="blue")
    #}#if
    return(p)
  }#if
  print(paste("A1HOM Background A2 threshold: ", a1hom_int))
  print(paste("A2HOM Background A1 threshold: ", a2hom_int))
}#function


 plotGTseqReads <- function(gtseq, locus){
  sub <- gtseq[gtseq$Locus==locus,]
  if (sum(sub$Total)>0){
    p <- ggplot(sub)+
      geom_histogram(aes(A1_reads+A2_reads))+
      theme_classic()+
      ggtitle(locus)+
      xlab("Reads per reaction")+
      ylab("Frequency")
    return(p)
  }#if
}#function



#par(mfrow=c(1,2))
for (l in unique(gtseq$Locus)){
  plotGTseqReads(gtseq, l)
  plotAlleleBalance(gtseq, l)
  #ggpubr::ggarrange(a,b)
}#for

