# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform differential expression analysis of each study, using limma
#' and available covariates: Gender, Age, HeartFailure, HTx (if different HF sampleshave different
#' status)
#' 
#' Inputs :  METAheart object with all experiments unified in a list with target and GEX info
#' 
#' Disclaimer: target files must include at least one column named HeartFailure

library(limma)

source("src/data_utils.R")
source("src/misc_utils.R")

METAheart = readRDS(file = "data/METAheart.rds") #main object

# Explore the covariates available:
# This will help to automatically define which studies require which type of correction

covariates_summary = lapply(METAheart, function(x){
   
  
  covariates = c("Age","Gender",
                  "HTx","DCM",
                  "Sample","HeartFailure",
                  "Disease","disease")
  
  red_targets = dplyr::filter(x$TARGETS,
                              HeartFailure == "yes")
   
   metadata = colnames(x$TARGETS)
   metadata_check = c("no","no","no","no","no")
   names(metadata_check) = c("Age","Gender",
                             "HTx", "DCM",
                             "Extra")
   
   if("Age" %in% metadata){
     metadata_check[1] = "yes"
   }
   
   if("Gender" %in% metadata){
     
     n_vars = length(unique(x$TARGETS$Gender))
     
     if(n_vars == 2){
       metadata_check[2] = "yes"
     }
     
   }
   
   if("HTx" %in% metadata){
     
     n_vars = length(unique(red_targets$HTx))
     
     if(n_vars == 2){
       metadata_check[3] = "yes"
     }
     
   }
   
   if("DCM" %in% metadata){
     
     n_vars = length(unique(red_targets$DCM))
     
     if(n_vars == 2){
       metadata_check[4] = "yes"
     }
     
   }
   
   if(sum(!(metadata %in% covariates))>0){
     metadata_check[5] = "yes"
   }
   
   return(metadata_check)
})

covariate_summary = data.frame(t(data.frame(covariates_summary)),
                               stringsAsFactors = F)
covariate_summary$ID = rownames(covariate_summary)


saveRDS(covariate_summary, file = "data/covariate_summary.rds")

# 1) Run special analysis with diabetes for GSE26887

covariate_summary = dplyr::filter(covariate_summary,
                                  !ID %in% "GSE26887")

HF_mat = METAheart$GSE26887$GEX
targets = METAheart$GSE26887$TARGETS

hf = factor(targets$HeartFailure, levels= c("no","yes"))
diabetes = factor(targets$Diabetes, levels= c("no","yes"))
gender = factor(targets$Gender, levels= c("female","male"))
age = targets$Age

design = model.matrix(~hf+diabetes+gender+age)
fit = lmFit(object = HF_mat, design = design)
fit2 = eBayes(fit)

HF_results_adj = as.data.frame(topTable(fit2,adjust.method = "BH",
                                        coef="hfyes",number = Inf)) %>% 
  rownames_to_column() %>%
  arrange(desc(abs(t))) %>% as_tibble()

colnames(HF_results_adj)[1] = "ID"

METAheart$GSE26887[["HF_limma"]] = HF_results_adj

# 2) Run special analysis for Translatome

covariate_summary = dplyr::filter(covariate_summary,
                                  !ID %in% "Translatome")

HF_mat = METAheart$Translatome$GEX
targets = METAheart$Translatome$TARGETS

hf = factor(targets$HeartFailure, levels= c("no","yes"))
htx = factor(targets$HTx, levels= c("yes","no"))
gender = factor(targets$Gender, levels= c("female","male"))
age = targets$Age

design = model.matrix(~hf+htx+gender+age)

#batch adjustment:

corfit = duplicateCorrelation(HF_mat,design,
                              block=targets$batch_sign)

fit = lmFit(HF_mat, design,
            block = targets$batch_sign,
            correlation = corfit$consensus.correlation)

fit2 = eBayes(fit)

HF_results_adj = as.data.frame(topTable(fit2,adjust.method = "BH",
                                        coef="hfyes",number = Inf)) %>% 
  rownames_to_column() %>%
  arrange(desc(abs(t))) %>% as_tibble()

colnames(HF_results_adj)[1] = "ID"

METAheart$Translatome[["HF_limma"]] = HF_results_adj

# 3) Run standard linear model wo/ additional covariates = 5 tests

ids = dplyr::filter(covariate_summary,
                    Age == "no" &
                    Gender == "no",
                    HTx == "no") %>%
      dplyr::select(ID)

covariate_summary = dplyr::filter(covariate_summary,
                                  !ID %in% ids[[1]])

METAheart[ids[[1]]] = lapply(METAheart[ids[[1]]], function(x){
  x[["HF_limma"]] = run_HFlimma(x$GEX, x$TARGETS)
  return(x)
})

# 4) Run the analysis to data sets that require control of gender only = 1 test

ids = dplyr::filter(covariate_summary,
                    Age == "no" &
                      Gender == "yes",
                    HTx == "no") %>%
  dplyr::select(ID)

covariate_summary = dplyr::filter(covariate_summary,
                                  !ID %in% ids[[1]])

METAheart[ids[[1]]] = lapply(METAheart[ids[[1]]], function(x){
  HF_mat = x$GEX
  targets = x$TARGETS
  
  hf = factor(targets$HeartFailure, levels= c("no","yes"))
  gender = factor(targets$Gender, levels= c("female","male"))
  
  design = model.matrix(~hf+gender)
  fit = lmFit(object = HF_mat, design = design)
  fit2 = eBayes(fit)
  
  HF_results_adj = as.data.frame(topTable(fit2,adjust.method = "BH",
                                          coef="hfyes",number = Inf)) %>% 
    rownames_to_column() %>%
    arrange(desc(abs(t))) %>% as_tibble()
  
  colnames(HF_results_adj)[1] = "ID"
  
  x[["HF_limma"]] = HF_results_adj
  
  return(x)
})


# 5) Run the analysis to data sets that require control of gender & age only = 6

ids = dplyr::filter(covariate_summary,
                    Age == "yes" &
                      Gender == "yes",
                    HTx == "no") %>%
  dplyr::select(ID)

covariate_summary = dplyr::filter(covariate_summary,
                                  ! ID %in% ids[[1]])

METAheart[ids[[1]]] = lapply(METAheart[ids[[1]]], function(x){
  HF_mat = x$GEX
  targets = x$TARGETS
  
  hf = factor(targets$HeartFailure, levels= c("no","yes"))
  gender = factor(targets$Gender, levels= c("female","male"))
  age = targets$Age
  
  design = model.matrix(~hf+gender+age)
  fit = lmFit(object = HF_mat, design = design)
  fit2 = eBayes(fit)
  
  HF_results_adj = as.data.frame(topTable(fit2,adjust.method = "BH",
                                          coef="hfyes",number = Inf)) %>% 
    rownames_to_column() %>%
    arrange(desc(abs(t))) %>% as_tibble()
  
  colnames(HF_results_adj)[1] = "ID"
  
  x[["HF_limma"]] = HF_results_adj
  
  return(x)
})

# 6) Run the analysis to data sets that require control of HTx only

ids = dplyr::filter(covariate_summary,
                    Age == "no" &
                      Gender == "no",
                    HTx == "yes") %>%
  dplyr::select(ID)

covariate_summary = dplyr::filter(covariate_summary,
                                  ! ID %in% ids[[1]])

METAheart[ids[[1]]] = lapply(METAheart[ids[[1]]], function(x){
  HF_mat = x$GEX
  targets = x$TARGETS
  
  hf = factor(targets$HeartFailure, levels= c("no","yes"))
  htx = factor(targets$HTx, levels= c("yes","no"))

  design = model.matrix(~hf+htx)
  fit = lmFit(object = HF_mat, design = design)
  fit2 = eBayes(fit)
  
  HF_results_adj = as.data.frame(topTable(fit2,adjust.method = "BH",
                                          coef="hfyes",number = Inf)) %>% 
    rownames_to_column() %>%
    arrange(desc(abs(t))) %>% as_tibble()
  
  colnames(HF_results_adj)[1] = "ID"
  
  x[["HF_limma"]] = HF_results_adj
  
  return(x)
})

# 7) Check if success

if(nrow(covariate_summary)==0) print("finished succesfully")

saveRDS(METAheart, file = "data/METAheart.rds")
