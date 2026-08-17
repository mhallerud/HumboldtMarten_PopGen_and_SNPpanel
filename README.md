# Humboldt_Marten_Genetic_Assessment_scripts

## Citation
Hallerud, M.A., K. Moriarty, K.L. Pilgrim, C.E. Eriksson, M.K. Schwartz, and T. Levi. 2026. Leveraging next-generation tools for genetic assessment and demographic monitoring in threatened and elusive Humboldt martens. *Evolutionary Applications* 19(7): e70277. [https://doi.org/10.1111/eva.70277](https://doi.org/10.1111/eva.70277).

## Contents:
### Data
* GBAS_Individual_Genotypes.csv: Consensus genotypes for each individual.
* GBAS_Raw_Genotypes.csv: Replicate-level genotypes for each sample.

### Code
* **1-SNP_discovery.sh**: ddRADseq SNP discovery with *de novo* Stacks
* **2-SNP_Filtering.sh**: SNP filtering for ddRADseq dataset.
* **2B-SNP_Filtering_GBAS_Panel.sh**: SNP filtering to identify candidate SNPs for use in GBAS panel.
* **3-PVCA_batch_effect.R**: Assessing batch effects via principal variance components analysis.
* **4-PCoAs_ddRADseq.R**: Running principal coordinates analysis (PCoAs) for ddRADseq data (various filtering stages).
* **5-ADMIXTURE.sh**: Running ADMIXTURE analysis.
* **5B-ADMIXTUREplots_and_evalAdmix.R**: Plotting ADMIXTURE and evalAdmix results.
* **6-fineRADstructure.sh**: Running fineRADstructure.
* **7-ddRADseq_PopGenMetrics.R**: Calculating population genetics metrics (e.g., Fst, Fis, Ho, Hs) and individual metrics (e.g., heterozygosity and inbreeding) from the ddRADseq dataset.
* **8-triangulaR.R**: Assessing putative hybridization with triangle plots.
* **9-NeEstimator.R**: Calculating effective population size with NeEstimator.
* (oops- no 10)
* **11-IQtree_phylogeny.sh**: Constructing phylogeny using IQtree-2.
* **12-mitogenome_nucleotide_diversity.sh**: Calculating nucleotide diversity from previously published mitogenomes.
* **13-GTseq_Genotyping.txt**: Processing raw GBAS sequencing data using the [GTseq-Pipeline](https://github.com/GTseq/GTseq-Pipeline/).
* **14-GBAS_SNP_Discovery.txt**: Discovering additional SNPs from within the raw GBAS sequencing data.
* **15-GBAS_pairwise_matches_individualID.R**: Calculating the number of pairwise matching SNPs, pairwise PIDsib, and pairwise # SNPs genotyped for GBAS sample genotypes.
* **16-GBAS_panel_assessment.R**: Plotting genotyping success rates and off-target rates, calculating PID metrics, etc.
* **17-RAD_GBAS_comparison.R**: Comparison of ddRADseq and GBAS inference for equivalent individuals.
* **18-Full_GBAS_Analyses.R**: Analyses based on the full GBAS dataset (validation + application).


## Additional Data Resources:
Minimally demultiplexed ddRADseq data is available in NCBI's Sequence Read Archive under [PRJNA1470952](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1470952/).
