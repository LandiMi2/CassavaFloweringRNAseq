# Cassava Flowering RNA-seq

Scripts for the RNA-seq analysis of three cassava genotypes (DSC120, DSC196, DSC296) sampled at three flowering time points (T0, T1, T2), with and without light supplementation (54 samples, 3 replicates each).

Paper: [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.08.06.743231v1)

## Scripts

| File | Description |
|---|---|
| `QCandMapping.sh` | Read quality control and mapping to the cassava genome (STAR) |
| `Design.R` | Loads the featureCounts matrix, filters low-count genes, builds the sample table, and runs DESeq2 (design `~ condition * geno * time`) |
| `PCA.R` | Principal component analysis of the samples |
| `DEGsGenoEffect.md` | Pairwise genotype contrasts and differentially expressed genes (adjusted p < 0.05, \|log2FoldChange\| > 1) |
| `DEGsLightEffect.md` | Light vs no-light contrasts and differentially expressed genes |
| `GO.md` | GO enrichment of the DEGs with clusterProfiler |

