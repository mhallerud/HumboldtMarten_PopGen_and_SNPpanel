#### CALCULATE MATCHING SNPS BETWEEN POPULATIONS ####
### PURPOSE: Assigning sample genotypes to individuals
library(reshape)

# function to calculate pairwise matches
CalcPairwiseMatches <- function(consensus, type="prop"){
  # set up df to hold outputs
  gensim <- data.frame(matrix(NA,nrow=length(unique(consensus$SampleID)),ncol=length(unique(consensus$SampleID))))
  rownames(gensim) <- colnames(gensim) <- unique(consensus$SampleID)
  # loop through each sample
  for (i in unique(consensus$SampleID)){
    #print(i)
    for (j in unique(consensus$SampleID)){
      #print(paste0(".....",j))
      # NOTE: ONLY USE COLUMNS WITH "MAIN" SNPS (i.e., NOT microhaplotypes)
      consensusub <- consensus[which(consensus$SampleID%in%c(i,j)), c(14:106)]
      if(i==j){
        consensusub <- rbind(consensusub,consensusub)
      }#if
      nasnps <- which(consensusub[1,]=="NA" | consensusub[2,]=="NA" | is.na(consensusub[1,]) | is.na(consensusub[2,]) | consensusub[1,]=="0" | consensusub[2,]=="0")
      if(length(nasnps)>0){
        consensusub <- consensusub[,-nasnps]
      }#if
      if(!is.data.frame(consensusub)){
        consensusub <- as.data.frame(matrix(consensusub, nrow=2, byrow=T))
      }#if
      sim <- sum(consensusub[1,]==consensusub[2,], na.rm=TRUE)
      if(type=="prop_match"){
        gensim[which(rownames(gensim)==i), which(colnames(gensim)==j)] <- sim / ncol(consensusub)
      }else if(type=="snp_match"){
        gensim[which(rownames(gensim)==i), which(colnames(gensim)==j)] <- sim
      }else if(type=="geno_count"){
        gensim[which(rownames(gensim)==i), which(colnames(gensim)==j)] <- ncol(consensusub)
      }else{
        print("'type' parameter must be one of 'prop_match', 'snp_match', or 'geno_count'")
        stop()
      }#ifelse
    }#for j
  }#for i
  return(gensim)
}#CalcPairwiseMatches

# this file has consensus genotypes for all samples that successfully genotyped
genos <- read.csv("GBAS_genotypes.csv")
consensus <- genos[which(genos$Filepath=="CONSENSUS"),]
dim(consensus)#480 consensus genotypes

# calculate pairwise proportion matching, # snps matching, and # snps genotyped in both samples
gensim <- CalcPairwiseMatches(consensus, type="snp_match")
gensim_prop <- CalcPairwiseMatches(consensus, type="prop_match")
gensim_count <- CalcPairwiseMatches(consensus, type="geno_count")

# plot distributions
hist(reshape::melt(gensim)[,2], main="Pairwise Genotype Similarity", xlab="# matching SNPs")
hist(reshape::melt(gensim_prop)[,2], main="Pairwise Genotype Similarity", xlab="% matching SNPs")
hist(reshape::melt(gensim)[,2], main="Pairwise Genotyping Success", xlab="# SNPs genotyped in both samples")

# save
write.csv(gensim, "MACA_genotypes_similarity.csv") # raw number of matching SNPs
write.csv(gensim_prop, "MACA_prop_gensim_.csv") # proportion of matching SNPs for SNPs genotyped in both samples
write.csv(gensim_count, "MACA_pairwise_geno.csv") # SNPs genotyped in both samples

### NOTE: The MACA_genotypes_similarity.csv file was next used to manually compare samples with high # matching SNPs