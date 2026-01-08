#!/bin/bash

# versions of software:
# PLINK v1.90b7.4 64-bit (18 Aug 2024)
# VCFtools (0.1.17)
# bcftools 1.21
# bwa 0.7.19-r1273
# seqkit 2.8.2; https://doi.org/10.1371/journal.pone.0163962; https://doi.org/10.1002/imt2.191


##### SNP FILTERING FOR REDUCED SNP PANEL ####
##### NOTE: THIS WAS CONDUCTED PRIOR TO RUNNING THE THIRD DDRADSEQ BATCH SO ONLY INCLUDED BATCHES 1-2 #####
##### INPUT: VCF FOLLOWS QUALITY FILTERING IN 2A-SNP_Filtering_Genomics, 
### populations.loci.fa is raw output from Stacks populations module with --fasta-loci flag



#### 0) INITIAL MISSINGNESS FILTERS
# iterative filtering to remove SNPs with high missingness
plink --file HighMissingnessRAD.vcf --geno 0.20 --out 11A_geno20 --recode -aec
plink --file 11A_geno20 --mind 0.40 --out 11B_geno20mind40 --recode -aec



#### 1) REMOVE SNPS IN FIRST 25 bp OF READ TO ALLOW FOR FWD PRIMER BINDING 
## NOTE: 25 bp is used here to allow for 20-bp of primer binding excluding the initial 5-bp RAD site
## NOTE: This didn't end up being a problem, but the full locus should be removed whenever a SNP occurs in the flanking region
## (this filtering allows loci through if there are additional SNPs outside of flanking regions)

## remove loci where SNP is within 25 bp of start
awk '{ if ($2>25) print $0}' 11C_geno30mind40geno20.vcf > 12A_QualityFiltered.snps.notFirst25.vcf
# NOTE: this file is missing its header so isn't a true vcf yet...

## fix CHROM field
# merge into single file
cat 12A_QualityFiltered.snps.notFirst25.vcf | tr ":" "\t" | awk 'FS="\t" {print $3}' > 12A_QualityFiltered.snps.notFirst25.chroms
paste 12A_QualityFiltered.snps.notFirst25.chroms 12A_QualityFiltered.snps.notFirst25.vcf | cut -f 1,3- > 12A_QualityFiltered.snps.notFirst25.fixed.vcf

# Fix VCF header before proceeding...
head 11C_geno30mind40geno20.vcf -n 7 > header.vcf
cat header.vcf 12A_QualityFiltered.snps.notFirst25.vcf > 12A_QualityFiltered.snps.notFirst25.fixed.vcf




#### 2) REMOVE LOCI WITH SNPs IN LAST 20 bp OF READ TO ALLOW FOR REV PRIMER BINDING ####
# make file with CLocus_ID in first column and sequence in second
grep '>' populations.loci.fa | cut -d_ -f2 > IDs.loci.fa
grep -v '>' populations.loci.fa > sequences.loci.fa

# calculate sequence lengths
awk '{print $1}' sequences.loci.fa | awk '{print length}' > loci.lengths
paste IDs.loci.fa loci.lengths sequences.loci.fa > lengths.loci.fa

# link loci lengths to VCFs, then filter for SNPs that aren’t within last 20 bp of the read
awk 'NR==FNR {h[$1] = $2; next} {print h[$1] "\t" $0}' lengths.loci.fa 12A_QualityFiltered.snps.notFirst25.fixed.vcf |
awk '{if ($1>$3+20) print $0}' | cut -f2- > 12B_QualityFiltered.snps.notLast20.vcf

# Fix VCF before proceeding...
cat header.vcf 12B_QualityFiltered.snps.notLast20.vcf > 12B_QualityFiltered.snps.notLast20.fixed.vcf

# make list of locus IDs for SNPs not in first 25 bp or last 20 bp
bcftools query -f '%ID\n' 12B_QualityFiltered.snps.notLast20.fixed.vcf > 12B_QualityFiltered.snps.notFlankingRegions.ids

# extract these...
vcftools --vcf 11C_geno30mind40geno20.vcf --snps 12B_QualityFiltered.snps.notFlankingRegions.ids --out 12B_QualityFiltered.snps.notFlankingRegions --recode




#### 3) REMOVE LOCI ALIGNING TO Y CHROMOSOME ####
## Align stacks loci to Y chromosome
# Index reference chrom
bwa index ../../../info/Mustela_erminea_Y_chrom.fna.gz
# Align
bwa mem ../../../../../info/Mustela_erminea_Y_chrom.fna.gz populations.loci.fa > populations.snps.Ychrom.sam


## Quality filtering of alignments
awk '{if($2!=4) print $0}' populations.snps.Ychrom.sam |
# remove any sequences with bitwise flag >=256 = 0x100 (these are secondary alignments or chimeric alignments, this is a standard threshold)
awk '{if ($2<256) print $0}' | 
# remove any sequences with MapQ score < 20 (this is roughly equivalent to requiring 99% certainty in the match)
awk '{if ($5>=20) print $0}' > populations.snps.Ychrom.filtered.sam


## Remove loci that align to Y
# Make list of IDs that align to Y
awk '{print $1}' populations.snps.Ychrom.filtered.sam | cut -d'_' -f 2 > populations.snps.Ychrom.filtered.loci

# Remove these from the loci file
seqkit grep -v -f populations.snps.Ychrom.filtered.loci populations.loci.fa > populations.noY.loci.fa 

# Extract SNPs in these loci from the raw VCF
## NOTE: These loci likely won't be in the filtered dataset since they're only in males and probably didn't meet missingess criteria!
awk 'NR==FNR {A[$1]=$3; next} $1 in A {print $0}' populations.snps.Ychrom.filtered.loci populations.snps.vcf > populations.snps.Yloci.vcf

# Remove these loci from the filtered VCF
awk '{print $3}' populations.snps.Yloci.vcf > populations.snps.Yloci.ids
vcftools --vcf 12B_QualityFiltered.snps.notFlankingRegions.recode.vcf --out 13A_Filtered_noYloci --exclude populations.snps.Yloci.ids --recode
#vcftools --vcf 12B_QualityFiltered.snps.notFlankingRegions.recode.vcf --out 13B_Filtered_Yloci --snps populations.snps.Yloci.ids --recode




#### 4) ALLELE FREQUENCY FILTERING FOR INDIVIDUAL ID ####
## filter for all SNPs with MAF>0.10 (within each coastal marten population)
## this should make the panel somewhat useful throughout Pacific marten range, 
## though we'll select from this subset to optimize specifically for coastal marten ID


## FILTERING FOR INITIAL PANEL (FIRST 50 SNPs)
## MAF ~0.20 ACROSS POPULATIONS
## THESE WERE BASED ONLY ON FIRST PLATE OF SAMPLES & ALSO PRIOR TO POP STRUCTURE ANALYSIS

## make list of locus IDs with MAF<0.20 in each pop, 0.1 for pops with low sample size
# skipping Lassen and Rockies because sample size is N=2 for each...
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --keep Siuslaw.txt --maf 0.20 --out 1A_snps.filtered.Siuslaw --recode
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --keep NCali.txt --maf 0.15 --out 1B_snps.filtered.NCali --recode
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --keep BlueMountains.txt --maf 0.20 --out 1C_snps.filtered.BlueMountains --recode
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --keep Cascades.txt --maf 0.10 --out 1D_snps.filtered.Cascades --recode
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --keep Trinidad.txt --maf 0.10 --out 1E_snps.filtered.Trinidad --recode

bcftools query -f '%ID\n' 1A_snps.filtered.Siuslaw.recode.vcf > 2A_Siuslaw.maf20.IDs
bcftools query -f '%ID\n' 1B_snps.filtered.NCali.recode.vcf > 2B_NCali.maf15.IDs
bcftools query -f '%ID\n' 1C_snps.filtered.BlueMountains.recode.vcf > 2C_BlueMountains.maf20.IDs
bcftools query -f '%ID\n' 1D_snps.filtered.Cascades.recode.vcf > 2D_Cascades.maf10.IDs
bcftools query -f '%ID\n' 1E_snps.filtered.Trinidad.recode.vcf > 2E_Trinidad.maf10.IDs

# combine MAF20 locus IDs into one list
cat 2A_Siuslaw.maf20.IDs 2B_NCali.maf15.IDs 2C_BlueMountains.maf20.IDs > 2_AllPops.maf20.ids

# tabulate # pops with MAF>0.20 for each locus ID, then extract locus IDs with MAF20 for all 4 pops
sort 2_AllPops.maf20.ids | uniq -c | sort -nr | awk '{if ($1==3) print $2}' > 2_AllPops.maf20.cleaned.ids

# extract these from VCF
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --snps 2_AllPops.maf20.cleaned.ids --out 3_snps.MAF20allpops --recode




### 4B) FILTERING FOR EXPANDED PANEL (+100 SNPs):
# filter for SNPs meeting MAF>0.10 across coastal Martens
vcftools --vcf 13A_Filtered_noYloci.recode.vcf --maf 0.10 --keep CoastalMartens.txt --out CoastalMartens.maf10 --recode 





#### 5) ALIGN TO MARTES FLAVIGULA GENOME ####
# subset stacks to given loci
grep -v "#" CoastalMartens.maf10.recode.vcf | awk '{print "CLocus_"$1}' > CoastalMartens.maf10.chroms
seqkit grep -v populations.loci.fa -f CoastalMartens.maf10.chroms > CoastalMartenMAF10.loci.fa


# Align Stacks loci to reference assembly
bwa mem ../info/Martes_flavigula_v1.fna.gz CoastalMartenMAF10.loci.fa -t 30 | samtools view -b | samtools sort > CoastalMartenMAF10_MartesFlavigula.bam

# Quality filtering of alignments
# remove any sequences with bitwise flag >=256 = 0x100 (these are secondary alignments or chimeric alignments, this is a standard threshold)
# remove any sequences with MapQ score < 20 (this is roughly equivalent to requiring 99% certainty in the match)
samtools view --excl-flags 0x800 -q 20 CoastalMartenMAF10_MartesFlavigula.bam | awk '{print $1}' > CoastalMartenMAF10_MartesFlavigula.filtered.loci

# extract these loci from the populations FA...
seqkit grep -v populations.loci.fa -f CoastalMartenMAF10_MartesFlavigula.filtered.loci | \
# convert multiline to single FA...
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' > CoastalMartenMAF10.MAFL.fa






##### 6) TRIM END OF READS #####
cutadapt CoastalMartenMAF10.MAMAaligned.loci.fa -a CCGAGATCGGAAG -o CoastalMartenMAF10.MAMA.trimmed.fa -j 10 




#### 7) CENSOR REPETITIVE ELEMENTS ####
# First, run the through Censor Repbase- save masked part as fasta: https://www.girinst.org/censor/
# select sequence source: Carnivora
# then, remove any reads with masked portions "X"
grep -v -B1 X CoastalMartenMAF10.CENSORmask.fa > CoastalMartenMAF10.CENSORed.fa

# finally, find these loci in VCF...
grep ">" CoastalMartenMAF10.CENSORed.fa | sed 's/>//g' > CoastalMartenMAF10.CENSORed.loci
grep -v "#" CoastalMartens.maf10.recode.vcf | awk '{print $1,$2,$3}' | cut -d":" > CoastalMartens.maf10.snps
awk -v RS='\r?\n' 'FNR==NR{arr[$0];next} ($1 in arr)' CoastalMartenMAF10.CENSORed.loci CoastalMartens.maf10.snps | awk '{print $3}' > CoastalMartenMAF10.CENSOR.ids
vcftools --vcf CoastalMartenMAF10.CENSORed.fa -snps CoastalMartenMAF10.CENSOR.ids --recode --out CoastalMartenMAF10.CENSOR




#### THIS VCF + POPULATIONS.LOCI.FA WAS THEN INPUT TO MULTIPLEX_WORMHOLE (https://github.com/mhallerud/multiplex_wormhole) FOR PRIMER DESIGN ####
## NOTE: per-population MAFs (NDunes, SDunes, NCali, SOreg) were further assessed in R using vcfR::maf function to assess designed primer sets prior to ordering!