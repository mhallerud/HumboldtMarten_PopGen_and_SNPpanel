# single commands were run on the CLQS cluster

# first, convert low-missingness VCF to PHYLIP using https://github.com/edgardomortiz/vcf2phylip
vcf2phylip.py -i LDprunedRAD.vcf -m 1 --output-prefix LDpruned

# use IQtree to ID and remove invariant sites... This will throw an error but print variant sites to LDpruned.min1.phy.varsites.phy
iqtree2 -s LDpruned.min1.phy -m MFP+ASC -B 1000

# step 1: use ModelFinder to identify the model with highest likelihood + build UNROOTED maximum-likelihood phylogeny
# -m: ModelFinder + ascertainment bias correction
# -B: 1000 ultra-fast bootstraps to quantify branching support
iqtree2 -s LDpruned.min1.phy.varsites.phy -m MFP+ASC -B 1000 --prefix LDpruned

# step 2: infer rooted tree based on DNA substitutions and quantify rootstrap support for different root positions
# -m: UNREST (i.e., unrestricted - most general non-reversible DNA substitution model)
# --root-test: re-roots existing tree on every branch and quantifies support
# -au: perform approximately-unbiased (AU) test for the tree
# -zb: 1000 rootstraps for quantifying root support and tests 
# -te: existing phylogeny
# --prefix: outfiles prefix 
iqtree2 -s LowMissingnessRAD.min1.phy.varsites.phy -m UNREST --root-test -au -zb 1000 -te LowMissingnessRAD.min1.phy.varsites.phy.treefile --prefix LowMissingnessRAD_UNRESTroottest

# consensus tree used in visualization = *rootstrap.contree. Bootstrap support in *roostrap.nex, test results in log file