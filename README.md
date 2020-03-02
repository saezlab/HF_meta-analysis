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
 
Raw data of each experiment used is not provided as it can be downloaded from their original publications. However, scripts showing how data was individually processed is provided [here](put_link)
 
***

### How to cite
> 

***

### Analyses & Scripts
#### Generation of list of data sets used in all analysis (already provided in Zenodo)
Script available [here](put_link).
ID modifications for manuscript [here](put_link).

#### Differential expression analysis
Analysis script available [here](put_link).

#### Differential expression analysis, gene level statistics visualization (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Generation of list of external data sets used (already provided in Zenodo)
Script available for external [here](put_link) and for fetal studies [here](put_link).

#### Data description: PCAs of all data sets and z-transformed data sets, t-SNE (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Data description: PCAs of gene-centered data (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Data description: DCM vs ICM for gene std matrices (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Replicability of individual studies and transfer learning classifier (Figure 2)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Robustness of replicability measurements (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Meta-analysis and generation of consensus ranking (Figure 3)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Selection of most informative genes from the consensus ranking
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Added value of meta-analysis
Analysis script available [here](put_link).

#### Gene level variability (Supplemental Figure #)
Analysis script available [here](put_link).
Figure script available [here](put_link).

#### Functional transcriptomics
Analysis script available [here](put_link).
Figure script available [here](put_link).



