#!/bin/bash
# tool versions:
# vcftools 0.1.17
# PLINK v1.90b5.2
# bcftools 1.21



### QUALITY FILTERING FOR PACIFIC MARTEN DDRADSEQ ###
#### NOTE: This pipeline largely follows the recommendations in O'leary et al. 2018 & McKinney et al. 2017



#### 1) RUN PER-SITE DEPTH CALCULATIONS ####
## pull per-site depth data from VCF
bcftools query -f '%CHROM %POS[\t  %DP]\n' populations.snps.vcf > bcftools.depths.log

## remove depths with NA (i.e., individuals where no SNP was present)
sed 's/\./NA/g' bcftools.depths.log > bcftools.depths.na.log

## take sums of depths per site - this should go up to column 42!
awk '{FS=OFS}{for(i=3;i<=99;i++) sum+=$i}{print $0, sum}' bcftools.depths.na.log > bcftools.depths.sum

## take counts of individuals where SNP was present (for calculating mean)
awk '{for(i=3;i<=99;i++) if($i!="NA") n[$1]+=1; print $0, n[$1]; n[$1]=0}' bcftools.depths.sum > bcftools.depths.n

## take mean of depths per site (sum/n from previous)
awk '{FS=OFS}{print $0, $100/$101}' bcftools.depths.n > bcftools.depths.mean

## take standard deviation of depths per site (based on mean and n previously calculated)
awk '{for(i=3;i<=99;i++) if($i!="NA") sumsq[$1]+=($i-$102)^2; if(($101-1)!=0) print $0, sqrt(sumsq[$1]/($101-1)); else print $0, Inf; sumsq[$1]=0}' bcftools.depths.mean > bcftools.depths.sd




#### 2) REMOVE LOW-QUALITY GENOTYPES ####
## remove loci with min depth <5 (per locus per individual) - can change minDP to 3 if using FreeBayes since uses info from whole pop to call SNPs...
# with 5 reads, the probability of falsely calling a homozygote is 6.25% (per locus) - for 4 reads: 12.5%, for 6 reads: 3.12%
vcftools --vcf populations.snps.subset.vcf --minDP 5 --out 1_populations.snps.minDP5 --recode

## remove loci with minor allele count <=2
vcftools --vcf 1_populations.snps.minDP5.recode.vcf --mac 2 --out 2_populations.snps.mac3 --recode

## mean depth (per site) > 8
# marten lowest 10th percentile=2.5x, 20th=6.7, 25th=8.1
vcftools --vcf 2_populations.snps.mac3.recode.vcf --min-meanDP 8 --out 3_populations.snps.meanDP8 --recode


## quality score > 20 (99% precision)
# extract SNPs with genotype quality (quality of SNP site) >=20
bcftools view -i 'GQ>=20' 3_populations.snps.meanDP8.recode.vcf > 4_populations.snps.Q20.vcf

# convert to plink for next step
vcftools --vcf 4_populations.snps.Q20.vcf --out 4_populations.snps.Q20 --plink

# fix plink file before moving on
mv 4_populations.snps.Q20.map 4_populations.snps.Q20.orig.map
awk '$1="un"' 4_populations.snps.Q20.orig.map > 4_populations.snps.Q20.map




#### 3) REMOVE MISSINGNESS #### 
## iteratively remove low quality loci and individuals - final params >70% genotyping rate and <50% missing per individual
plink --file 4_populations.snps.Q20 --geno 0.50 --out 5A_populations.snps.geno50 --recode --allow-extra-chr
plink --file 5A_populations.snps.geno50 --mind 0.80 --out 5B_populations.snps.mind80 --recode --allow-extra-chr
plink --file 5B_populations.snps.mind80 --geno 0.40 --out 5C_populations.snps.geno30 --recode --allow-extra-chr
plink --file 5C_populations.snps.geno30 --mind 0.50 --out 5D_populations.snps.mind50 --recode vcf --allow-extra-chr



#### 4) MITIGATE ALLELE DROPOUT / COVERAGE EFFECTS #####
## remove loci with high variance in read depth among individuals - this could indicate allele dropout
# Marten Threshold: 99th percentile of depth SDs = 37.5
#rm MaxSD_DP37.IDs -f
#awk '{if(($NF)>37.5) print $1 ":" $2}' bcftools.depths.stats > MaxSD_DP37.IDs
#vcftools --vcf 5D_populations.snps.mind50.vcf --exclude MaxSD_DP37.IDs --out 6_populations.snps.maxSD_DP37 --recode
## SKIPPED THIS STEP BECAUSE 0 SNPS REMOVED

## NOTE: FIX SAMPLE IDs IN THIS OUTPUT BEFORE PROCEEDING!!!



#### 5) REMOVE LOCI OUTSIDE OF HWE IN MOST POPULATIONS ####
## NOTE: Populations only considered here if they had reasonable sample sizes. Trinidad + NCali merged.
## before running this, make keep files with sample IDs per population
## each output VCF will have SNPs per population
vcftools --vcf 5C_mind50.vcf --keep NDunes.txt --out 6A_populations.snps.NDunes --plink
vcftools --vcf 5C_mind50.vcf --keep SDunes.txt --out 6B_populations.snps.SDunes --plink
vcftools --vcf 5C_mind50.vcf --keep Lassen.txt --out 6C_populations.snps.Lassen --plink
vcftools --vcf 5C_mind50.vcf --keep BlueMountains.txt --out 6D_populations.snps.BlueMountains --plink
vcftools --vcf 5C_mind50.vcf --keep SouthernOregon.txt --out 6E_populations.snps.SouthernOregon --plink
vcftools --vcf 5C_mind50.vcf --keep NCali.txt --out 6F_populations.snps.NCali --plink
vcftools --vcf 5C_mind50.vcf --keep Denver.txt --out 6G_populations.snps.Denver --plink

## calculate HWE p-values per site per pop
plink --file 6A_populations.snps.NDunes --hardy --out 7A_populations.snps.NDunes
plink --file 6B_populations.snps.SDunes --hardy --out 7A_populations.snps.SDunes
plink --file 6C_populations.snps.Lassen --hardy --out 7A_populations.snps.Lassen
plink --file 6D_populations.snps.BlueMountains --hardy --out 7A_populations.snps.BlueMountains
plink --file 6E_populations.snps.SouthernOregon --hardy --out 7A_populations.snps.SouthernOregon
plink --file 6F_populations.snps.NCali --hardy --out 7A_populations.snps.NCali
plink --file 6G_populations.snps.Denver --hardy --out 7A_populations.snps.Denver

## calculate Benjamini-Hochberg (1995) p-value correction for multiple tests
# Benjamini-Hochberg procedure:
# 1. Put the individual p-values in ascending order.
# 2. Assign ranks to the p-values. For example, the smallest has a rank of 1, the second smallest has a rank of 2.
# 3. Calculate each individual p-value’s Benjamini-Hochberg critical value, using the formula (i/m)Q, where:
#       i = the individual p-value’s rank,
#       m = total number of tests,
#       Q = the false discovery rate (a percentage, chosen by you).
# 4. Compare your original p-values to the critical B-H from Step 3; find the largest p value that is smaller than the critical value.
        # assign # rows to variable
        # sort by HWE p-value
        # add rank, taking into account equivalent values
        # calculate BH p-value based on rank
        # export IDs for corrected p-vals less than alpha

## For each population, the output file will hold SNP IDs that are NOT in HWE based on B-H procedure
for pop in NDunes SDunes NCali SouthernOregon BlueMountains Lassen Denver
do
  	ROWS=$(wc -l 7A_populations.snps.${pop}.hwe | awk '{print $1}')
        sort 7A_populations.snps.${pop}.hwe -k 9 | \
        awk 'BEGIN{i=0}{if($9>prev) i++; print $0, i; prev=$9}' | \
        awk -v ROWS="$ROWS" '{print $0, ($10/ROWS)*0.01}' | \
        awk '{if ($9<0.05) print $2}' > ${pop}.HWD.IDs
done

## Following Pearman et al. 2022, will only exclude loci out of HWE in all populations.
## Due to low sample sizes in some populations, we'll only require HWD in 5 (of 7) pops.
## HWE alpha value: 0.01 - conservative due to small sample size in many pops

# count number of populations that each SNP is not in HWE
cat *.HWD.IDs > 7_allpopsHWD.IDs
sort 7_allpopsHWD.IDs | uniq -c | sort -nr | awk '{if ($1>=3) print $2}' > 7_HWD3pops.IDs

# remove SNPs not in HWE in 5 or more populations
vcftools --vcf 5C_mind50.vcf --exclude 7_HWD3pops.IDs --out 7_populations.snps.hwe3pops --recode



#### 6) REMOVE PUTATIVE PARALOGS ####
## remove super high  read depths: mean + 2*SD = 45 
# NOTE: stacks also does this!
vcftools --vcf 7_populations.snps.hwe3pops.recode.vcf --max-meanDP 45 --out 8_populations.snps.maxDP45 --plink

# fix plink file
mv 8_populations.snps.maxDP45.map 8_populations.snps.maxDP45.orig.map
awk '$1="un"' 8_populations.snps.maxDP45.orig.map > 8_populations.snps.maxDP45.map


## remove excess heterozygosity (ObsHet > 0.65 when HWE>0.05)
# calculate allele freqs & HWE p-value
plink --file 8_populations.snps.maxDP45 --hardy --out 8_populations.snps.maxDP45 --allow-extra-chr

# calculate Benjamini-Cochberg correction on p-values
# HWE alpha value: 0.05
awk 'BEGIN{i=1}{print $0, i}END{i+=1}' 8_populations.snps.maxDP45.hwe | \
awk '{print $0, ($10/NR)*0.05}' > 8_populations.snps.maxDP45.BC.hwe

# make list of loci with Obs_Het>0.60 [column 6] --only if not in HWE! [column 9]
# HWE alpha value 0.05 - again to be conservative
awk '{if($7>0.60 && $9>0.05) print $2}' 8_populations.snps.maxDP45.hwe > Het60.IDs

# extract loci
plink --file 8_populations.snps.maxDP45 --exclude Het60.IDs --out 9_populations.snps.Het60 --recode vcf --allow-extra-chr



#### 7) REMOVE BIASED SNPS #####
## remove loci with large deviations from allele balance: <=0.2 or >=0.8 homozygosity:heterozygosity
# extract allelic depth info from original SNP file
bcftools query -f '%CHROM %POS [%AD\t]\n' populations.snps.vcf > AlleleDepths.log

# calculate allele balance (major allele depth / total site depth)
sed 's/\./0,0/g' AlleleDepths.log | \
        sed 's/ /\t/g' | \
        tr ',' '\t' | \
        awk '{for (i=1; i<=16; i++) if($(2*i+1)!=0 || $(2*i+2)!=0) print $1 ":" $2, $(2*i+1) / ($(2*i+2)+$(2*i+1)); else print $1 ":" $2, "."}' > AlleleBalance.log

# make list of sites with AB>0.80 or AB<0.20 - these are likely have allele dropout!
awk '{if($2!=1 && $2!=0 && $2!="." && ($2<0.20 || $2>0.80)) print $1}' AlleleBalance.log > UnbalancedAlleles.IDs

# remove SNPs
vcftools --vcf 9_populations.snps.Het60.vcf --exclude UnbalancedAlleles.IDs --out HighMissingnessRAD_withDups --recode
mv HighMissingnessRAD_withDups.recode.vcf HighMissingnessRAD_withDups.vcf #rename


# convert to plink & fix plink file
vcftools --vcf 10_AB_duplicates --out 10_AB_duplicates --plink
mv 10_AB_duplicates.map 10_AB_duplicates.orig.map
awk '$1="un"' 10_AB_duplicates.orig.map > 10_AB_duplicates.map



#### 7B) REMOVE DUPLICATES ####
## Duplicates.txt: sampleIDs of replicates
vcftools --vcf HighMissingnessRAD.recode.vcf --remove Duplicates.txt --out HighMissingnessRAD --recode
mv HighMissingnessRAD.recode.vcf HighMissingnessRAD.vcf #rename

#### THIS COMPRISES THE "HIGH-MISSINGNESS" DATASET USED IN ANALYSES ####



#### 8) CREATE LOW-MISSINGNESS DATASET #####
# ultimately, missingness doesn't have a huge effect on patterns observed in data-- aim for genotyping / missingness per individual rate >50%
# consider --geno filtering per population
plink --vcf HighMissingnessRAD.vcf --const-fid --geno 0.05 --out 11A_geno05 --recode --allow-extra-chr
plink --file 11A_geno05 --mind 0.10 --mac 1 --out LowMissingnessRAD --recode --allow-extra-chr



#### 9) LINKAGE DISEQUILIBRIUM PRUNING ####
# do this last to avoid losing SNPs due to errors!
# --indep-pairwise (window size in SNPs, size of window shift at each step, Pearson's correlation threshold)
plink --vcf HighMissingnessRAD.vcf --out 11_LD50560 --indep-pairwise 50 5 0.6 --allow-extra-chr
vcftools --vcf HighMissingnessRAD.vcf --snps 11_LD50560.prune.in --out LDprunedRAD --recode


# also create low missingness pruned dataset...
plink --vcf LowMissingnessRAD.vcf --out 12_geno05mind10LD50560 --indep-pairwise 50 5 0.6 --allow-extra-chr
vcftools --vcf LowMissingnessRAD.vcf --snps 12_geno05mind10LD50560.prune.in --out 12_geno05mind10LD50560 --recode
