#!/bin/bash
####------------RUN ADMIXTURE MODELS FOR DDRADSEQ------------####
VCF=$1 #load input VCF as first argument


### INPUT VCFs RUN:
## A) full LD-pruned dataset: LDprunedRAD.vcf
## B) LD-pruned data with first-relatives removed: 11B
# to remove close relatives: 
plink2 --pedmap x --make-king-table -aec
awk '{if($8>0.2) print $0}' plink2.kin0  
# choose one per pair to remove...
## C) one sample per geographic region
## D) A but only Humboldt Martens


# extract basename of input VCF
NAME="$(basename $VCF | cut -d'.' -f1)"

# convert VCF to BED file
vcftools --vcf ${VCF} --plink --out ${NAME}
mv ${NAME}.map ${NAME}.orig.map 
awk '$1="0"' ${NAME}.orig.map > ${NAME}.map
plink --file ${NAME} --out ${NAME} --make-bed -aec


# run ADMIXTURE for K 1-12
for K in 1 2 3 4 5 6 7 8 9 10 11 12
do
	# set up folder
  	[ -d K${K} ] || mkdir K${K}
        cd K${K}

	# run admixture
	admixture --cv=10 ../${NAME}.bed $K > "${NAME}"_"${K}".log
        cd ..

	# run evalAdmix
	evalAdmix/evalAdmix -plink ${NAME} -fname K${K}/${NAME}.${K}.P -qname K${K}/${NAME}.${K}.Q -o evalAdmix_${NAME}_k${K}.out
done

# print cross-validation error for all K-values run
grep 'CV error' K*/"${NAME}"*.log
