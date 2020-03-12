# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Description: In this script we model gene specific
# variability with 2-way ANOVAs

source("HGEX_src/data_utils.R") #general functions 
source("HGEX_src/misc_utils.R")
METAheart = readRDS(file = "HGEX_data/METAheart.rds") #main object

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

marrays = c("GSE57345","GSE42955",
            "GSE1869","GSE3585",
            "GSE26887","GSE5406")

load("./HGEX_data/dictionaryIDs.ro")
new_ids = dictionary %>% dplyr::filter(GEO_ID %in% marrays)
marrays = new_ids$newID

meta_targets = meta_targets %>% mutate(Tech = ifelse(ExpID %in% marrays,
                                                     "microarray","rnaseq"))

# Modelling each factor

fisher_rank = run_fisher_meta(meta_list = METAheart ,
                              n_missing = 4)

genes = names(fisher_rank)

# Performing a 2-way anova with interaction term (HF/study)

anova_study = run_anovastats(numeric_matrix = meta_gex[genes,],
                             targets = meta_targets) #For sup

anova_study_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                   targets = meta_targets) #For sup



# Performing a 2-way anova with interaction term (HF/technology)

anova_tech = run_anovastats(numeric_matrix = meta_gex[genes,],
                            targets = meta_targets,
                            factor_b = "Tech")

anova_tech_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                  targets = meta_targets,
                                  factor_b = "Tech") 

# Performing a 2-way anova with interaction term (HF/gender)

gender_targets = filter(meta_targets, !is.na(Gender))

anova_gender = run_anovastats(numeric_matrix = meta_gex[genes,gender_targets$grl_id],
                              targets = gender_targets,
                              factor_b = "Gender")

anova_gender_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,gender_targets$grl_id],
                                    targets = gender_targets,
                                    factor_b = "Gender") 

# Performing a 2-way anova with interaction term (HF/Age)

age_breaks = c(0,35,65,100)

age_targets =  meta_targets %>% filter(!is.na(Age)) %>% 
  mutate(Age = .bincode(x = Age,
                        breaks = age_breaks,
                        include.lowest = TRUE))

anova_age = run_anovastats(numeric_matrix = meta_gex[genes,age_targets$grl_id],
                           targets = age_targets,
                           factor_b = "Age")

anova_age_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,
                                                                 age_targets$grl_id],
                                 targets = age_targets,
                                 factor_b = "Age")

# Performing a 2-way anova with interaction term (HF/Age)

anova_explanted = run_anovastats(numeric_matrix = meta_gex[genes,],
                             targets = meta_targets,
                            factor_b = "HTx") 

anova_explanted_scale = run_anovastats(numeric_matrix = meta_gex_scale[genes,],
                                       targets = meta_targets,
                                       factor_b = "HTx")

# Supplementary tables

WriteXLS(x = c("anova_study",
               "anova_tech",
               "anova_gender",
               "anova_age",
               "anova_explanted"), 
         ExcelFileName = "HGEX_data/paper_sup/confounding_anova_raw.xlsx",
         SheetNames = c("Study",
                        "Tech",
                        "Gender",
                        "Age",
                        "HTx")
)


WriteXLS(x = c("anova_study_scale",
               "anova_tech_scale",
               "anova_gender_scale",
               "anova_age_scale",
               "anova_explanted_scale"), 
         ExcelFileName = "HGEX_data/paper_sup/confounding_anova.xlsx",
         SheetNames = c("Study",
                        "Tech",
                        "Gender",
                        "Age",
                        "HTx")
)


# Figure processing

anova_study = mutate(anova_study, 
                     data_type = "raw") %>% 
  filter(stats=="etasq") 

anova_study_scale = mutate(anova_study_scale, 
                           data_type = "gene_centered") %>%
  filter(stats=="etasq")

anova_study_res = bind_rows(anova_study,anova_study_scale)

####
anova_tech = anova_tech %>% mutate(data_type = "raw") %>% 
  filter(stats=="etasq")

anova_tech_scale = anova_tech_scale %>% mutate(data_type = "gene_centered") %>% 
  filter(stats=="etasq")

anova_tech_res = bind_rows(anova_tech,anova_tech_scale)

####

anova_gender = anova_gender %>% mutate(data_type = "raw") %>% 
  filter(stats=="etasq")

anova_gender_scale = anova_gender_scale %>% mutate(data_type = "gene_centered") %>% 
  filter(stats=="etasq")

anova_gender_res = bind_rows(anova_gender,
                             anova_gender_scale)

####
anova_age = anova_age %>% mutate(data_type = "raw") %>% 
  filter(stats=="etasq")

anova_age_scale = anova_age_scale %>% mutate(data_type = "gene_centered") %>% 
  filter(stats=="etasq")

anova_age_res = bind_rows(anova_age,anova_age_scale)

####
anova_explanted = anova_explanted %>% mutate(data_type = "raw") %>% 
  filter(stats=="etasq")

anova_explanted_scale = anova_explanted_scale %>% mutate(data_type = "gene_centered") %>% 
  filter(stats=="etasq")

anova_explanted_res  = bind_rows(anova_explanted,
                                 anova_explanted_scale)

# Plotting objects

anova_confounding = list("study" = anova_study_res,
                         "tech" = anova_tech_res,
                         "gender" = anova_gender_res,
                         "age" = anova_age_res,
                         "HTx" = anova_explanted_res)

saveRDS(anova_confounding, 
        file = "HGEX_data/figure_objects/confounding_anova.rds")

































