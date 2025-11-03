Repository Overview

Horlings et al. dataset. A preprocessed version of the Horlings array-CGH dataset is included in this repository so users can directly apply the TAaCGH algorithm and the extended pipeline described in this study. The dataset is redistributed in processed form for reproducibility and educational use only.

TCGA dataset. Open-access files were obtained from the NCI Genomic Data Commons (GDC). In the respective folder, we provide scripts to download the official files directly from the GDC API and to reproduce all preprocessing into the matrices used by this study. We include only small, de-identified example subsets in this repository to keep the repo light and avoid duplicating large third-party datasets.

METABRIC dataset. Clinical and copy-number data were obtained from the Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016) study available at the cBioPortal for Cancer Genomics. The files brca_metabric_clinical_data.tsv and METABRICCNADiscovery.txt were used in conjunction to extract and analyze copy-number alterations in Chromosome 5 corresponding to the PURPL (LINC01021) locus. These data are redistributed here only in minimal, de-identified subset form for demonstration purposes; users should download the full datasets directly from the cBioPortal source.

SCAN-B (GSE96058) dataset. Normalized RNA-seq expression data and corresponding clinical metadata were obtained from the NCBI Gene Expression Omnibus (GEO; accession GSE96058). The scripts in this repository demonstrate how to retrieve, clean, and process the official data files following the R data package workflow described in https://12379monty.github.io/GSE96058/. To keep this repository lightweight and in compliance with data redistribution policies, only small, de-identified example subsets are included. The full dataset must be downloaded directly from the original GEO source. Mutation data were incorporated using files downloaded from the SCAN-B Mutation Explorer: https://oncogenomics.bmc.lu.se/MutationExplorer/.



File Structure and Functionality

Files ending in _survival.R perform Kaplan–Meier survival analyses using either copy number or gene expression data. Adjusted p-values are calculated using the Benjamini–Hochberg (BH) false discovery rate correction. Files ending in _data_download.R are intended to retrieve official datasets from their original repositories (such as GEO, Mutation Explorer, or GDC) and/or prepare them for use in this study. The Data/(Horlings;METABRIC;SCAN-B;TCGA) folder contains de-identified example datasets that demonstrate how to run the provided scripts without redistributing protected or large-scale data.



Main Analysis: Extended TAaCGH Pipeline

The central script extending the TAaCGH workflow is cnv_region_probe2gene.R. After identifying significant regions using TAaCGH, this script locates genes within those regions that are statistically significant. Running this pipeline generates a Results/ directory organized by breast cancer subtype or characteristic (e.g., Luminal A, Luminal B, Basal, HER2, PI3K, etc.).

Each subtype folder includes a table of BH-adjusted (FDR) p-values for each probe in the region. Optional Bonferroni-corrected values (FWER) are also provided but are not emphasized in this study. Each folder also contains a file of patient profiles relevant to that subtype and a list of significant probes from the previously identified TAaCGH regions, annotated by BH or Bonferroni significance.

For each probe identified as significant by BH correction, a dedicated subfolder is created. Within each probe-specific folder, adjacent probes are used to define the genomic interval surrounding the probe. The Ensembl gene browser is then queried to retrieve annotation data—including HGNC symbols, Ensembl gene IDs, gene biotypes, and genomic coordinates—for all genes located within or near the interval defined by the adjacent probes. The resulting files summarize the genes found in each region and indicate which genes lie exactly at the probe location. Each folder also contains boxplots illustrating copy-number or expression differences by subtype (Luminal A, Luminal B, Basal, HER2, PI3K, etc.), along with a companion file reporting the statistical results for these plots. In addition, a biotype summary file describes the composition and characteristics of genes located between adjacent probes.

