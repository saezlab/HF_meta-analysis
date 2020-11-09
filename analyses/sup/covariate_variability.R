# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform PCAs of individual data.
#' Then we assign a percentage of variable for each useful predictor variable
#' 
#' 
#' 
#' 

#Inputs and dependencies
source("src/data_utils.R") #general functions 
source("src/misc_utils.R")

#Covariate summary
covariate_summary = readRDS("data/covariate_summary.rds")


# Dictionary of studies
METAheart = readRDS(file = "data/METAheart.rds") #main object

# First - scale data + PCA
METAheart = lapply(METAheart, function(x){
  x[["GEX_norm"]] = PLIER::rowNorm(x[["GEX"]])
  pca_res = summary(prcomp(t(x[["GEX_norm"]]),
                       center = T,scale. = T))
  x[["PCA"]] = list("PCs" = pca_res$x,
                    "prop_var" = pca_res$importance[2,])
  return(x)
})

# Second - calculate the prop_variance of each covariate

tidyPCA = function(x){
  targets = x[["TARGETS"]]
  
  pcs = as.data.frame(x[["PCA"]][["PCs"]]) %>%
    rownames_to_column(.,"Sample") %>%
    tidyr::pivot_longer(-Sample,
                 names_to = "PC") %>%
    dplyr::left_join(targets,by = "Sample") %>%
    dplyr::arrange(PC) %>%
    dplyr::group_by(PC) %>%
    nest()
  
  return(pcs)
}

get_propvar_all = function(data,
                       pos_tests = c("Age","Gender","HeartFailure","Diabetes")){
  
  useful_columns = set_names(colnames(data)[colnames(data) %in% pos_tests])
  
  covariate_tests = lapply(useful_columns, function(covariate){
    
    model_data = data[,c("value",covariate)]
    colnames(model_data) = c("value","covariate")
    
    if(is.numeric(model_data$covariate) |
       length(unique(model_data$covariate)) > 1){
      
    results_lm = summary(lm(value~covariate,model_data))
    
    if(results_lm$coefficients[2,4]<0.05){
      return(TRUE)
    }else{
      return(FALSE)
    }
    
    }else{
      return(NA)
    }
  }) %>% enframe() %>% unnest()
  
  return(covariate_tests)
  
}

study_variance = lapply(METAheart, function(x){
  
  pcs = tidyPCA(x) %>% 
    mutate(assoc_test = map(data,get_propvar_all)) %>% 
    dplyr::select(PC,assoc_test) %>%
    unnest() 
  
  pcs_T = pcs %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::filter(value == TRUE)
  
  pcs_F = pcs %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::filter(value == FALSE)
  
  pcs_N = pcs %>%
    dplyr::filter(is.na(value))
  
  if(nrow(pcs_T) > 0){
    
    tidy_prop_var = tibble("PC" = names(x$PCA$prop_var), 
                           "prop_var" = x$PCA$prop_var)
    
    pcs_T = left_join(pcs_T, tidy_prop_var) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(name) %>%
      summarise(total_var = sum(prop_var))
    
    covars_w_data = unique(pcs_T$name)
    
    pcs_all = pcs_T
    
    #Complete data:
    
    #Covariates that were tested but didn't associate with PCs = 0
    
    pcs_F = dplyr::filter(pcs_F, !name %in% covars_w_data)
    
    if(nrow(pcs_F)>0){
      pcs_F_vals = tibble("name" = unique(pcs_F$name), 
                          "total_var" = 0)
      pcs_all = bind_rows(pcs_all, pcs_F_vals)
    }
    
    #Covariates that are reported but have no more than 1 class = NULL
    
    pcs_N = dplyr::filter(pcs_N, !name %in% covars_w_data)
    
    if(nrow(pcs_N)>0){
      pcs_N_vals = tibble("name" = unique(pcs_N$name), 
                          "total_var" = NA)
      pcs_all = bind_rows(pcs_all, pcs_N_vals)
    }
    #Covariates that weren't tested = NULL
    
    not_tested = c("Age","Gender","HeartFailure","Diabetes")
    not_tested = not_tested[!not_tested %in% unique(pcs$name)]
    if(length(not_tested) > 0 ){
      pcs_absent_vals = tibble("name" = not_tested, 
                               "total_var" = NA)
      pcs_all = bind_rows(pcs_all, pcs_absent_vals)
    }
    
    return(pcs_all)
  }else{
    
    #Covariates that were tested but didn't associate with PCs = 0
    covars_w_data = unique(pcs_F$name)
    pcs_all = tibble("name" = unique(pcs_F$name), 
                        "total_var" = 0)
    
    #Covariates that are reported but have no more than 1 class = NULL
    pcs_N = dplyr::filter(pcs_N, 
                          !name %in% covars_w_data)
    
    if(nrow(pcs_N)>0){
      pcs_N_vals = tibble("name" = unique(pcs_N$name), 
                          "total_var" = NA)
      pcs_all = bind_rows(pcs_all, pcs_N_vals)
    }
    
    #Covariates that weren't tested = NULL
    
    not_tested = c("Age","Gender","HeartFailure","Diabetes")
    not_tested = not_tested[!not_tested %in% unique(pcs$name)]
    if(length(not_tested) > 0 ){
      pcs_absent_vals = tibble("name" = not_tested, 
                               "total_var" = NA)
      pcs_all = bind_rows(pcs_all, pcs_absent_vals)
    }
    
    return(pcs_all)
  }
}) %>% enframe("study") %>% unnest()


## Plot ## To be moved depending on the result

experiment_size = sort(unlist(lapply(METAheart,
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

original_features = ggplot(study_variance,aes(y = factor(study,
                                     levels = names(experiment_size)),
                          x = name,
                          fill = total_var)) + 
  geom_tile() + theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 90, hjust=1)) +
  ylab("Studies") + 
  xlab("Covariates") +
  scale_fill_gradient(
    low = "black",
    high = "yellow",
    limits = c(0,0.6))


pdf("data/figures/sup/covariate_variability.pdf", height = 8, width = 6)

plot(original_features)

dev.off()


saveRDS(study_variance,
        file = "data/figure_objects/covariate_var_all.rds")






















