## A Meta-Analysis on the Transcriptomic Landscape of end-stage Heart Failure

### Abstract

**Aims:** 
Transcriptomic studies have contributed to fundamental knowledge of myocardial remodeling in human heart failure (HF). However, the agreement on the crucial genes in HF is limited and systematic efforts to integrate evidences of multiple patient cohorts are lacking.  Here we aimed to provide an unbiased consensus transcriptional signature of human end-stage HF by comprehensive comparison and analysis of publicly available datasets. 

**Methods and Results:** 
We curated and uniformly processed 16 public transcriptomic studies of left ventricular samples from 263 healthy and 653 failing human hearts. Transfer learning approaches revealed conserved disease patterns across all studies independent of technical differences. We meta-analyzed the dysregulation of 14041 genes to extract a consensus signature of HF. Estimation of the activities of 343 transcription factors, 14 signalling pathways, and 182 micro RNAs, as well as the enrichment of 5998 biological processes confirmed the established aspects of the functional landscape of the disease and revealed novel ones. We provide all results in a free public resource [RefHF](https://saezlab.shinyapps.io/hgex_app/) to facilitate further use and interpretation of the results. We exemplify usage by deciphering fetal gene reprogramming and tracing myocardial origin of the plasma proteome biomarkers in HF patients.

**Conclusion:** 
We demonstrated the feasibility of combining transcriptional studies from different HF patient cohorts. In our compendium we provide a robust and consistent collection of molecular markers of end-stage HF that may guide the identification of novel targets with diagnostic or therapeutic relevance.

<img src="man/figures/graphical_abstract.png" align="center" width="800">

***

### Availabilty of data
The datasets supporting the conclusions of this publication are available at Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3564179.svg)](here_put_link)

From Zenodo you can download these (zipped) folders: 

 * `data` - Contains the directory with processed data sets used for the analysis and to generate the figures
 
Please deposit the unzipped folders in the root directory of this R-project.
 
 **Exceptions:**
 
Raw data of each experiment used is not provided as it can be downloaded from their original publications. However, scripts showing how data was individually processed is provided [here](https://github.com/saezlab/HF_meta-analysis/tree/master/data_processing/scripts)
 
***

### How to cite
> 

***

### Analyses & Scripts

#### Generation of list of data sets used in all analysis (already provided in Zenodo)
Script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main_objects/make_metaheart.R).

#### General description of studies (Figure 1)
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/sample_info_size.R).

#### Gene coverage (Supplemental Figure 2)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/gene_coverage.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gene_coverage_figs.R).

#### Differential expression analysis
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/de_analysis.R).

#### Differential expression analysis, gene level statistics visualization (Supplemental Figure 3)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/deg_stats.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/deg_stats.R).

#### Generation of list of external data sets used (already provided in Zenodo)
Script available for external [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main_objects/make_external_metaheart.R) and for fetal studies [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main_objects/make_fetal_metaheart.R).

#### Test of marker genes in all experiments (Supplemental Figure 4)
Script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/HF_marker_genes.R).

#### Data description: PCAs of all data sets and z-transformed data sets, t-SNE (Supplemental Figure 6)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/general_variability.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gen_var_figs.R).

#### Data description: PCAs of gene-centered data (Supplemental Figure 7)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/gene_centered_analysis.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gcentered_figs.R).

#### Data description: DCM vs ICM for gene std matrices (Supplemental Figure 8)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/dcm_vs_icm.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/dcm_vs_icm_figs.R).

#### Replicability of individual studies and transfer learning classifier (Figure 2)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/study_comparison.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/reproducibility_figs.R).

#### Robustness of replicability measurements (Supplemental Figure 9)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/robustness_glist_size.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/robustness_es_ds.R).

#### Meta-analysis and generation of consensus ranking (Figure 3)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/get_metaranking.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/meta_main.R).

#### Gradient of information in the meta-analysis (Supplemental Figure 10)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/genes_best_performance.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/best_perf_figs.R).

#### Added value of meta-analysis
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/added_value.R).

#### Gene level variability (Supplemental Figure 11)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/gene_variability.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gene_variability_anova.R).

#### Extrapolation of disease score (Supplemental Figure 12)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/ds_fetal_external_studies.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/ds_fetal_external_studies_plot.R).

#### Functional transcriptomics (Figure 4)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/functional_analysis.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/funcomics_tiles.R).

#### Exploration of Consensus Signature (Figure 5, Supplemental Figure 13)
Proteomic analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/validation_proteomic.R).
Fetal response analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/validation_fetal.R).
Main Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/validation_plotting.R).
Supplemental Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/supp_validation.R).

