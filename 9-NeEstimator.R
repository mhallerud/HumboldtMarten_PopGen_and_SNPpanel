#devtools::install_github("https://github.com/zakrobinson/RLDNe")
library(RLDNe)
library(adegenet)
#library(hierfstat)

# read in data
vcf <- read.vcfR("LowMissingnessRAD.vcf")
#popmap <- read.csv("../../SampleLocations/Marten_Hets_Fhom_coords.csv", stringsAsFactors=FALSE, fileEncoding='Latin1', header=T)

# grab population IDs 
pop(genind) <- unlist(lapply(indNames(genind), function(x) strsplit(x, "_")[[1]][1]))
unique(pop(genind))
# merge border with N Cali
pop(genind)[pop(genind=="Border")] <- "NCali"

# function to set up NEestimator files
run_NEstimator_LD <- function(prefix, vcf=vcf, pops=pop(genind)){
  # load genepop
  # NOTE: These were created by using the populations module from Stacks: populations -V <.vcf> -O . --genepop
  vcf_sub <- vcf[,c(1,which(pops%in%prefix)+1)]# subset to population
  vcf_sub <- vcf_sub[which(maf(vcf_sub)[,4]>0),]# remove fixed alleles
  vcf_sub
  coma <- vcfR2genind(vcf_sub, return.alleles=FALSE)
  pop(coma) <- rep(prefix, length(indNames(coma)))
  
  # convert to format
  genos <- as.data.frame(tab(coma))
  colnames(genos) <- str_replace_all(colnames(genos), fixed(".1"), "_1")
  colnames(genos) <- str_replace_all(colnames(genos), fixed(".0"), "_2")
  #genos <-genos[,-5]
  genotypes <- alleles2genotypes(df=genos, allele_cols=1:dim(genos)[2], allelesAsIntegers=TRUE)
  genotypes[genotypes==11] <- 12 # 1-1 (1 minor allele, 1 major allele)
  genotypes[genotypes==20] <- 22 # 2-0 (both alleles major)
  genotypes[genotypes==2] <- 11 # 0-2 (both alleles minor)
  
  gp_file <- write_genepop_zlr(loci=genotypes, 
                               pops=pop(coma),
                               ind.ids=indNames(coma),
                               folder="",
                               filename =paste0(prefix,"_output.txt"),
                               missingVal=NA, 
                               ncode=2,
                               diploid=T)
  
  param_files <- NeV2_LDNe_create(input_file=gp_file$Output_File,
                                  param_file=paste0(prefix,"_Ne_params.txt"),
                                  NE_out_file=paste0(prefix,"_Ne_out.txt"))

  run_LDNe(LDNe_params=param_files$param_file)
}#function run_NEstimator_LD



run_NEstimator_LD("NDunes")
run_NEstimator_LD("SDunes")
run_NEstimator_LD("SOreg")
run_NEstimator_LD("NCali")
run_NEstimator_LD("Lssn")
run_NEstimator_LD("BMT")
run_NEstimator_LD("CO")
