#!/bin/bash

### SNP Discovery Code for Pacific Marten ddRADseq ###
## using Stacks v2.65
## each step was run with multi-threading on OSU's CQLS cluster


## 1) demultiplex raw sequencing data
# set up folder structure
mkdir 01_Demultiplexed; cd 01_Demultiplexed
mkdir pstI_mspI_plate1
mkdir pstI_mspI_plate2
mkdir pstI_mspI_plate3

# run demultiplexing separately for each plate  
process_radtags -1 /nfs2/hts/nextseq/221116_VH00571_202_AACGYJYM5/L1/lane1-s001-index--GBS0334_S1_R1_001.fastq.gz \
-2 /nfs2/hts/nextseq/221116_VH00571_202_AACGYJYM5/L1/lane1-s001-index--GBS0334_S1_R2_001.fastq.gz \
-P -o pstI_mspI_plate1 -b ../info/barcodes_Plate1.txt -r -q -w 0.10 -c --score-limit 10 --renz-1 pstI --renz-2 mspI

process_radtags -1 /nfs2/hts/nextseq/230523_VH00571_273_AACTFLFM5/L1/lane1-s001-index--GBS0339_S1_R1_001.fastq.gz \
-2  /nfs2/hts/nextseq/230523_VH00571_273_AACTFLFM5/L1/lane1-s001-index--GBS0339_S1_R2_001.fastq.gz \
-P -o pstI_mspI_plate2 -b info/barcodes_Plate2.txt -r -q -w 0.15 -c --score-limit 10 --renz-1 pstI --renz-2 mspI

process_radtags -1 /nfs2/hts/nextseq/241203_VH00571_500_AAGGWGKM5/L1-GBS-1mm/lane1-s001-GBS0350_S1_R1_001.fastq.gz \
-2 /nfs2/hts/nextseq/241203_VH00571_500_AAGGWGKM5/L1-GBS-1mm/lane1-s001-GBS0350_S1_R2_001.fastq.gz \
-P -o pstI_mspI_plate3 -b ../info/barcodes_Plate3.txt -r -q -w 0.15 -c --score-limit 10 --renz-1 pstI --renz-2 mspI



## 2) kmer filtering for overly abundant and overly rare kmers
for s in `ls -1 pstI_mspI_plate*/*fq.gz | cut -d '.' -f1 | uniq`
do
	dir="$(dirname $s)"
	[ -d ${dir}/Kmer_Filtered ] || mkdir ${dir}/Kmer_Filtered
	kmer_filter -1 ${s}.1.fq.gz -2 ${s}.2.fq.gz -i 'gzfastq' -o ${dir}/Kmer_Filtered --rare --abundant
done



## 3) de novo Stacks pipeline
# set up folder structure
cd ..
#mkdir 02_Parameter_Optimization
mkdir 03_DeNovo_Stacks/; cd 03_DeNovo_Stacks
mkdir 02_Stacks/Marten/ -p
mkdir 03_SNP_Calls/Marten -p
  

## 3A) build loci from each sample for single-end reads
## NOTE: -m and -M parameters were selected based on parameter optimization using samples in plate1
id=1
for i in `ls -1 ../01_Demultiplexed_FQ/pstI_mspI_plate*/Marten/Kmer_Filtered/*.fq.gz | awk '{print substr($0, 1, length($0)-8 )}' | uniq`
do
  	ustacks -f "$i".1.fq.gz -o 02_Stacks/Marten -i $id -m 4 -M 2 -p 10
        let "id+=1"
done

# remove ".1" in filenames
for i in `ls 02_Stacks/Marten -1`
do
       out=$(ls 02_Stacks/Marten/$i | sed s/'\.1'//)
       #echo $out
       mv 02_Stacks/Marten/$i "$out"
done


## 3B) Build catalog of variable loci across samples
# NOTE: before this stp, move any samples that failed QC (based on CQLS) into a separate folder so that these aren't included in the catalog
mkdir 02_Stacks/Marten/FailedQC 
for i in MA-16-BMT8 Sbeach_M6A M14B Gbeach_M6A
do
       mv 02_Stacks/Marten/$i* 02_Stacks/Marten/FailedQC
done


# NOTE: Also move duplicates to a separate folder to avoid duplication in the catalog
mkdir 02_Stacks/Marten/Duplicates 
for i in ODFW_COMA_6856 COMA_6856_2 COMA_6856_3 COMA_6860_2 COMA_6860_3 MACA31_2
do
       mv 02_Stacks/Marten/$i* 02_Stacks/Marten/Duplicates
done

# build catalog of loci available in the population from samples in the pop map
## NOTE: -n parameter was set based on parameter optimization using samples from plate1
## NOTE: COMBINED_marten_popmap_catalog excludes duplicates and failed samples vs. the catalog used below
cstacks -P 02_Stacks/Marten -M ../info/COMBINED_marten_popmap_catalog.txt -n 1 -p 60 --report-mmatches


## 3C) Match samples against the catalog and call SNPs
# NOTE: Move FailedQC and Duplicate samples back into main folder before running this
mv 02_Stacks/Marten/FailedQC/* 02_Stacks/Marten
mv 02_Stacks/Marten/Duplicates/* 02_Stacks/Marten

# match samples in the popmap against the catalog
sstacks -P 02_Stacks/Marten -M ../info/COMBINED_marten_popmap_dunesplit.txt -p 60

# transpose data to store by locus and add in paired end reads
# NOTE: raw reads for plates 2-3 were copied into pstI_mspI_plate1/Marten/Kmer_Filtered to avoid having to run this separately
tsv2bam -P 02_Stacks/Marten -M ../info/COMBINED_marten_popmap_dunesplit.txt --pe-reads-dir ../01_Demultiplexed_FQ/pstI_mspI_plate1/Marten/Kmer_Filtered -t 60

# build PE contig from pop data, align reads per sample, and call variant sites in the population
gstacks -P 02_Stacks/Marten -M ../info/COMBINED_marten_popmap_dunesplit.txt -t 60 -O 03_SNP_Calls/Marten --write-alignments --gt-alpha 0.05 --details


## 3D) Run populations to make VCFs and other outputs
# run populations on known populations, subspecies, and sex
mkdir 03_SNP_Calls/Marten/populations
mkdir 03_SNP_Calls/Marten/sex
mkdir 03_SNP_Calls/Marten/subspecies

populations -P 03_SNP_Calls/Marten -M ../info/COMBINED_marten_popmap_dunesplit.txt -O 03_SNP_Calls/Marten/populations --vcf --plink --hwe --fstats --fasta-loci --genepop --radpainter --fasta-samples
populations -P 03_SNP_Calls/Marten -M ../info/COMBINED_marten_popmap_sex.txt -O 03_SNP_Calls/Marten/sex --vcf --plink --hwe --fstats --fasta-loci
populations -P 03_SNP_Calls/Marten -M ../info/COMBINED_marten_popmap_ssp.txt -O 03_SNP_Calls/Marten/subspecies --vcf --plink --hwe --fstats --fasta-loci
