# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de
# Description: In this script we model gene specific
# variability with 2-way ANOVAs

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

library(WriteXLS)


METAheart = lapply(METAheart, function(x){
  x[["GEX_norm"]] = PLIER::rowNorm(x[["GEX"]])
  return(x)
})

# Generating a unified data set

meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample,HeartFailure,Gender,Age,ExpID,HTx) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_"))

meta_gex = get_complete_gex(meta_list = METAheart,
                            complete_targets = meta_targets,
                            gex_key = "GEX")

meta_gex_scale = get_complete_gex(meta_list = METAheart,
                                  complete_targets = meta_targets,
                                  gex_key = "GEX_norm")

# Annotating RNAseq and Microarray

marrays = c("GSE76701","GSE57345","GSE42955",
            "GSE1869","GSE3585","GSE26887",
            "GSE5406","GSE16499")

load("./data/dictionaryIDs.ro")
new_ids = dictionary %>% dplyr::filter(GEO_ID %in% marrays)
marrays = new_ids$newID

meta_targets = meta_targets %>% mutate(Tech = ifelse(ExpID %in% marrays,
                                                     "microarray","rnaseq"))

# Modelling each factor

fisher_rank = run_fisher_meta(meta_list = METAheart ,
                              n_missing = length(METAheart) - 10)

genes = names(fisher_rank)

# Performing a 2-way anova with interaction term (HF/study)

anova_study_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                   targets = meta_targets) #For sup



# Performing a 2-way anova with interaction term (HF/technology)

anova_tech_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                  targets = meta_targets,
                                  factor_b = "Tech") 

# Performing a 2-way anova with interaction term (HF/gender)

gender_targets = filter(meta_targets, !is.na(Gender))

anova_gender_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,gender_targets$grl_id],
                                    targets = gender_targets,
                                    factor_b = "Gender") 

# Performing a 2-way anova with interaction term (HF/Age)

age_breaks = c(0,35,65,100)

age_targets =  meta_targets %>% filter(!is.na(Age)) %>% 
  mutate(Age = .bincode(x = Age,
                        breaks = age_breaks,
                        include.lowest = TRUE))

anova_age_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,
                                                                 age_targets$grl_id],
                                 targets = age_targets,
                                 factor_b = "Age")

# Performing a 2-way anova with interaction term (HF/Age)

anova_explanted_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                       targets = meta_targets,
                                       factor_b = "HTx")

WriteXLS(x = c("anova_study_scale",
               "anova_tech_scale",
               "anova_gender_scale",
               "anova_age_scale",
               "anova_explanted_scale"), 
         ExcelFileName = "data/paper_sup/confounding_anova.xlsx",
         SheetNames = c("Study",
                        "Tech",
                        "Gender",
                        "Age",
                        "HTx")
)


anova_confounding = list("study" = anova_study_scale,
                         "tech" = anova_tech_scale,
                         "gender" = anova_gender_scale,
                         "age" = anova_age_scale,
                         "HTx" = anova_explanted_scale)

anova_confounding = map(anova_confounding, function(x){
  
  mutate(x,data_type = "gene_centered") %>%
    dplyr::filter(stats=="etasq")
  
} )

saveRDS(anova_confounding, 
        file = "data/figure_objects/confounding_anova.rds")

































