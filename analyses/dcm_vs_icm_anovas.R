# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de
# Description: In this script we model gene specific
# variability with 2-way ANOVAs for ICM vs DCM comparison

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
library(cowplot)

METAheart = readRDS(file = "data/METAheart.rds") #main object

library(WriteXLS)

fisher_rank = run_fisher_meta(meta_list = METAheart ,
                              n_missing = length(METAheart) - 10)

genes = names(fisher_rank)

METAheart = lapply(METAheart, function(x){
  x[["GEX_norm"]] = PLIER::rowNorm(x[["GEX"]])
  return(x)
})

# Generating a unified data set

meta_targets = get_tibble_union(METAheart,"TARGETS") %>% 
  dplyr::select(Sample,HeartFailure,Gender,Age,ExpID,HTx,DCM) %>% 
  mutate(grl_id = paste(Sample,ExpID,sep = "_"))

meta_gex = get_complete_gex(meta_list = METAheart,
                            complete_targets = meta_targets,
                            gex_key = "GEX")

meta_gex_scale = get_complete_gex(meta_list = METAheart,
                                  complete_targets = meta_targets,
                                  gex_key = "GEX_norm")

# Modelling each factor

run_anovastatsDCM = function(numeric_matrix, targets, 
         factor_a = "HeartFailure", factor_b = "DCM"){
  
  library(sjstats)
  library(fastmatch)
  
  anova_res = apply(numeric_matrix, 1, function(x, targets){
    
    gene_i = x
    factor_a_vect = factor(targets[[factor_a]])
    factor_b_vect = factor(targets[[factor_b]])
    gene_aov = aov(gene_i ~ factor_a_vect + factor_b_vect)
    aov_stats = sjstats::anova_stats(gene_aov) 
    
  },targets = targets)  %>% bind_rows(.id = "gene") %>% as_tibble() %>% 
    group_by(gene) %>% gather(stats,value,-(gene:term))  %>% spread(term,value) %>%
    ungroup()
  
  colnames(anova_res)[fmatch(c("factor_a_vect","factor_b_vect"),
                             colnames(anova_res))] = c(factor_a,factor_b)
  
  return(anova_res %>% dplyr::filter(stats == "etasq") %>%
           pivot_longer(c(HeartFailure,DCM),names_to = "covariate") %>%
           dplyr::select(-Residuals))
}

# First all the data

allDCM_res = run_anovastatsDCM(numeric_matrix = meta_gex[genes,],targets = meta_targets)

allDCM_res %>% dplyr::filter(stats == "etasq") %>%
  pivot_longer(c(HeartFailure,DCM),names_to = "covariate") %>%
  dplyr::select(-Residuals) %>%
  dplyr::filter(gene %in%  genes[1:500]) %>%
  ggplot(aes(x= covariate, y = value)) + geom_boxplot() + theme_classic()

  dplyr::mutate(top500 = ifelse(gene %in%  names(genes)[1:500],"yes","no"))
  
# Then each individual data set separately
  
all_data = map(METAheart, function(x){
  pivot_wider(enframe(x))
}) %>% enframe(name = "study") %>% unnest()  %>%
  group_by(study) %>%
  dplyr::mutate(useful = map(TARGETS, function(x){
    
    if("DCM" %in% colnames(x)){
      
      hf_targets = dplyr::filter(x, HeartFailure == "yes")
      
      if(sum(table(hf_targets$DCM) >= 3) == 2){
        
        return(TRUE)
      } else{
        return(FALSE)
      }
      
    }else{return(FALSE)}
  })) %>% unnest(useful) %>%
  dplyr::filter(useful == TRUE) %>%
  dplyr::mutate(GEX = map(GEX, function(x){
    
    useful_genes = genes[genes %in% rownames(x)]
    return(x[useful_genes,])
    })) %>%
  dplyr::mutate(aov_results = map2(GEX,TARGETS,run_anovastatsDCM)) %>%
  ungroup()

all_aov = all_data %>%
  dplyr::select(study,aov_results)

saveRDS(all_aov, file = "./data/figure_objects/dcm_aov_all.rds")

# Perform wilcox test on useful genes

wilcox_results = all_aov %>%
  group_by(study) %>%
  mutate(wilcox = map(aov_results,function(df){
    
    filtered_df = df %>%
      dplyr::filter(gene %in%  genes[1:4000]) %>%
      dplyr::select(-stats) %>% 
      pivot_wider(names_from = covariate,values_from = value)
    
    stest = wilcox.test(filtered_df$HeartFailure,
                        filtered_df$DCM,
                        alternative = "greater",
                        paired = T) %>%
      broom::tidy()
    
    return(stest)
  })) %>% unnest(wilcox) %>%
  dplyr::mutate(plot_object = pmap(list(study,aov_results,p.value), function(x,y,z){
    
    filtered_df =  y %>%
      dplyr::filter(gene %in%  genes[1:4000]) %>%
      dplyr::mutate(covariate = ifelse(covariate == "DCM",
                                       "etiology","HF"))
    
    
    plt_object = ggplot(filtered_df,aes(x= covariate, y = value)) + geom_boxplot() + theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 8)) +
      ggtitle(paste0(x, " Wilcox p-value = ", round(z,2))) +
      ylab("")
    
    return(plt_object)
  }))
  
AOV_DCM_plots = plot_grid(wilcox_results$plot_object[[1]] +
            ylab("prop. of explained variance"),
          wilcox_results$plot_object[[2]],
          wilcox_results$plot_object[[3]],wilcox_results$plot_object[[4]],
          wilcox_results$plot_object[[5]] +
            ylab("prop. of explained variance"),
          wilcox_results$plot_object[[6]],
          wilcox_results$plot_object[[7]],wilcox_results$plot_object[[8]],
          ncol = 4)

pdf(file = "data/figures/sup/aov_dcm.pdf", width = 9, height = 5)
plot(AOV_DCM_plots)
dev.off()

