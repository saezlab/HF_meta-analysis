## A Meta-Analysis on the Transcriptomic Landscape of end-stage Heart Failure

### Abstract

**Background**

During heart failure (HF) the myocardium undergoes ventricular remodeling, which is accompanied by changes in gene expression. Transcriptomic studies have contributed to fundamental knowledge of this remodeling. However, many studies do not agree on the differentially expressed marker genes. We hypothesized that a comprehensive meta-analysis of available studies would help to address this incongruence and provide a more unified biological understanding of the transcriptome in HF. 

**Methods**

We curated and uniformly processed 14 public transcriptomic data sets consisting of left ventricular samples from 265 healthy and 642 failing hearts. We assessed their comparability by applying classifiers based on transfer learning. We then combined all datasets to extract a consensus gene signature. Combining this consensus signature with biological prior knowledge, we estimated transcription factor and signalling pathway activities.

**Results**

Although single studies reported highly dissimilar marker genes, transfer knowledge based approaches revealed that disease patterns are conserved across studies. We derived a gene ranking that reflects consistent molecular hallmarks of HF. The ranking is strongly enriched in fibrosis-related gene sets and contains significant footprints of numerous transcription factors including Nanog, Sox2, Pbx3, Mef2; pathways including Jak-Stat and miRNAs including MIR-206, MIR-514. Results are provided at https://saezlab.shinyapps.io/hgex_app/.

**Discussion**

We demonstrated the feasibility of combining transcriptional studies from different technologies, years and centers. We extracted a consensus gene signature as a valuable resource to understand common molecular patterns in HF. To our knowledge, this report represents the largest meta-analysis of transcriptional HF studies to date.

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

#### Differential expression analysis
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/de_analysis.R).

#### Differential expression analysis, gene level statistics visualization (Supplemental Figure 2)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/deg_stats.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/deg_stats.R).

#### Generation of list of external data sets used (already provided in Zenodo)
Script available for external [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main_objects/make_external_metaheart.R) and for fetal studies [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main_objects/make_fetal_metaheart.R).

#### Test of marker genes in all experiments (Supplemental Figure 3)
Script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/HF_marker_genes.R).

#### Data description: PCAs of all data sets and z-transformed data sets, t-SNE (Supplemental Figure 5)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/general_variability.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gen_var_figs.R).

#### Data description: PCAs of gene-centered data (Supplemental Figure 6)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/gene_centered_analysis.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gcentered_figs.R).

#### Data description: DCM vs ICM for gene std matrices (Supplemental Figure 7)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/dcm_vs_icm.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/dcm_vs_icm_figs.R).

#### Replicability of individual studies and transfer learning classifier (Figure 2)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/study_comparison.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/reproducibility_figs.R).

#### Robustness of replicability measurements (Supplemental Figure 8)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/robustness_glist_size.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/robustness_es_ds.R).

#### Meta-analysis and generation of consensus ranking (Figure 3)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/get_metaranking.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/meta_main.R).

#### Gradient of information in the meta-analysis (Supplemental Figure 9)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/sup/genes_best_performance.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/best_perf_figs.R).

#### Added value of meta-analysis
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/added_value.R).

#### Gene level variability (Supplemental Figure 10)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/gene_variability.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/sup/gene_variability_anova.R).

#### Functional transcriptomics (Figure 4)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/functional_analysis.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/funcomics_tiles.R).

#### Extrapolation of disease score (Figure 5)
Analysis script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/main/DS_Fetal_External.R).
Figure script available [here](https://github.com/saezlab/HF_meta-analysis/blob/master/analyses/figures/main/DS_Fetal_External_plotting.R).


