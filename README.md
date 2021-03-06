# Code and supplemental data for MOSTWAS: Multi-Omic Strategies for Transcriptome-Wide Association Studies (Bhattacharya *et al*, 2020)

1. Supplemental Data: Models and results from the manuscript.
  - ROS/MAP contains (1) a folder FinalModels_1MB with all final MOSTWAS ROS/MAP models, (2) AD_IAGP_reuslts.tsv with TWAS results for IGAP GWAS summary statistics for risk of Alzheimer's disease, (3) MDD_PGC_results.tsv with TWAS results for major depressive disorder risk using PGC summary statistics, and (4) MDD_UKBBGWAS_results.tsv with TWAS results for major depressive disorder risk using UKBB GWAX summary statistics
  - TCGA-BRCA contains (1) a folder FinalModels_1MB with all final MOSTWAS TCGA-BRCA models and (2) BRCAsurivval_ICOGs_results_update.tsv with TWAS results for ICOGs GWAS summary statistics for breast cancer-specific survival
  - Simulation contains (1) power_size.tsv with simulation results for power and size of distal-added last test and (2) sims_power_r2.tsv with simulation results for power and model performance of MOSTWAS models
  - Comparison With BGW-TWAS contains (1) compareROSMAP_inPSYCHencode.tsv with results of external predictive performance of MOSTWAS and BGW-TWAS models trained in ROS/MAP and imputed into PSYCHencode data and (2) compareTCGA.tsv with results of external predicitive performance of MOSTWAS and BGW-TWAS models trained in TCGA-BRCA and imputed into a hold-out set from TCGA-BRCA

2. Code used in manuscript.
  - SimulationPipeline.R: simulation from TCGA
  - Sample_MeTWASCode.R: sample code to run MeTWAS
  - Sample_DePMACode.R: sample code to run DePMA
  - Sample_BGW-TWAS_Code.sh: sample code used for BGW-TWAS in TCGA-BRCA data
