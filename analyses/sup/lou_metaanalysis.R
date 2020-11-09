# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we make a LOU crossvalidation
#' of the HF-CS

source("src/data_utils.R") #general functions 
source("src/misc_utils.R")
METAheart = readRDS(file = "data/METAheart.rds") #main object

# 0. Processing object

# 1. Original Meta-analysis

meta_analysis = run_fisher_meta(meta_list = METAheart,
                                n_missing = length(METAheart) - 10)

experiments = names(METAheart)
names(experiments) = experiments
experiment_size = sort(unlist(lapply(METAheart, 
                                     function(x) ncol(x$GEX))),
                       decreasing = T)

# 1. What if we take one study out from the whole collection

experiments = names(METAheart)
names(experiments) = experiments

LOU_rankings = get_LOU_metarankings(experiments = experiments,
                                    METAheart = METAheart,
                                    n_missing = length(METAheart) - 10)

LOU_rankings = process_LOU(LOU_rankings = LOU_rankings)

jaccard_LOU_rankings = get_jaccard_LOU(LOU_rankings_tibble = LOU_rankings,genelist = "top1000")

order_study = c(names(experiment_size), "none")

ggplot(jaccard_LOU_rankings, 
       aes(x = factor(StudyA,
                      levels = order_study), 
           y = factor(StudyB,
                      levels = order_study), 
           fill = JaccardIx)) + geom_tile() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


jaccard_LOU_rankings %>%
  dplyr::filter( StudyB == "none",
                ! StudyB == StudyA) %>%
  dplyr::mutate(top4 = ifelse(StudyA %in% names(experiment_size[1:4]),
                              "Largest4","Else")) %>%
  wilcox.test(JaccardIx~top4,.)
  
jaccard_LOU_rankings %>%
  dplyr::filter( StudyB == "none",
                 ! StudyB == StudyA) %>%
  dplyr::mutate(top4 = ifelse(StudyA %in% names(experiment_size[1:4]),
                              "Largest4","Else"),
                study_label = ifelse(StudyA %in% names(experiment_size[1:4]),
                                     StudyA,"")) %>%
  ggplot(aes(x = top4, y = JaccardIx,
             label = study_label)) + geom_point() +
  theme_classic() + ggrepel::geom_text_repel(size =3)


all_jaccard = jaccard_LOU_rankings %>%
  dplyr::filter( StudyB == "none",
                 ! StudyB == StudyA) %>%
  dplyr::mutate(top4 = ifelse(StudyA %in% names(experiment_size[1:4]),
                              "Largest4","Else"),
                study_label = ifelse(StudyA %in% names(experiment_size[1:4]),
                                     StudyA,""))

all_jaccard %>% group_by(top4) %>% summarise(mean(JaccardIx))

saveRDS(jaccard_LOU_rankings,
        file = "data/figure_objects/jaccard_LOU_rankings.rds")


# 2. What happens if we use only RNAseq?

# What happens if we use only RNAseq
# Get marrays IDs
marrays = c("GSE76701","GSE57345","GSE42955",
            "GSE1869","GSE3585","GSE26887",
            "GSE5406","GSE16499")

load("./data/dictionaryIDs.ro")
dictionary = dictionary %>% 
  dplyr::mutate(tech = ifelse(GEO_ID %in% marrays,
                              "marray","rnaseq"))

experiments = dictionary %>% dplyr::filter(tech == "rnaseq") %>% pull(newID)
names(experiments) = experiments

LOU_rankings_rnaseq = get_LOU_metarankings(experiments = experiments,
                                    METAheart = METAheart,
                                    n_missing = 2)

LOU_rankings_rnaseq = process_LOU(LOU_rankings = LOU_rankings_rnaseq)

jaccard_LOU_rankings_rnaseq = get_jaccard_LOU(LOU_rankings_tibble = LOU_rankings_rnaseq)

order_study = c(names(experiment_size[experiments]), "none")

ggplot(jaccard_LOU_rankings_rnaseq, 
       aes(x = factor(StudyA,
                      levels = order_study), 
           y = factor(StudyB,
                      levels = order_study), 
           fill = JaccardIx)) + geom_tile() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

# Only microarray

experiments = dictionary %>% dplyr::filter(tech != "rnaseq") %>% pull(newID)
names(experiments) = experiments

LOU_rankings_rnaseq = get_LOU_metarankings(experiments = experiments,
                                           METAheart = METAheart,
                                           n_missing = 2)

LOU_rankings_rnaseq = process_LOU(LOU_rankings = LOU_rankings_rnaseq)

jaccard_LOU_rankings_rnaseq = get_jaccard_LOU(LOU_rankings_tibble = LOU_rankings_rnaseq,genelist = "top1000")

order_study = c(names(experiment_size[experiments]), "none")

ggplot(jaccard_LOU_rankings_rnaseq, 
       aes(x = factor(StudyA,
                      levels = order_study), 
           y = factor(StudyB,
                      levels = order_study), 
           fill = JaccardIx)) + geom_tile() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


#Without big microarray

# What happens if we take all rnaseq away

get_LOU_metarankings = function(experiments, METAheart, n_missing){
  
  lou_meta_ranking = map(experiments, function(experiment){
    
    studies_used = experiments[experiment != experiments]
    
    meta_analysis = run_fisher_meta(meta_list = METAheart[studies_used],
                                    n_missing = n_missing)
    
    results = tibble("gene" = names(meta_analysis),
                     "BH_p_value" = meta_analysis,
                     "ranking" = 1:length(meta_analysis),
                     "relative_ranking" = 1 - (1:length(meta_analysis)/length(meta_analysis)),
                     "n_genes" = length(meta_analysis))
    
    return(results)
    
  })
  
  lou_meta_ranking[["none"]] =tibble("gene" = names(meta_analysis),
                                 "BH_p_value" = meta_analysis,
                                 "ranking" = 1:length(meta_analysis),
                                 "relative_ranking" = 1 - (1:length(meta_analysis)/length(meta_analysis)),
                                 "n_genes" = length(meta_analysis))
 
  return(lou_meta_ranking)
   
}

process_LOU = function(LOU_rankings){
  LOU_rankings_tibble = enframe(LOU_rankings,name = "excluded",value = "meta_analysis") %>%
    mutate(sign_genes = map(meta_analysis, function(x){
      gene_list = dplyr::filter(x,BH_p_value <= 0.05) %>%
        pull(gene)
      return(gene_list)
    })) %>%
    mutate(top1000 = map(meta_analysis, function(x){
      gene_list = dplyr::slice(x,1:1000) %>%
        pull(gene)
      return(gene_list)
    })) %>%
    mutate(ngenes = map(meta_analysis, function(x){
      return(nrow(x))
    })) %>% unnest(ngenes)
  
  return(LOU_rankings_tibble)
}


get_jaccard_LOU = function(LOU_rankings_tibble, genelist = "sign_genes"){
  
  experiments = set_names(LOU_rankings_tibble$excluded)
  study_genelist = set_names(LOU_rankings_tibble[[genelist]], experiments)
  
  # Make a matrix/data-frame with pairwaise comparisons
  jaccard_res_sign =  enframe(lapply(experiments, function(set_a){
    genes_a = study_genelist[[set_a]]
    j_ix_a = lapply(experiments, function(set_b){
      
      genes_b = study_genelist[[set_b]]
      
      #Jaccard Index
      j_ix = length(intersect(genes_a,genes_b))/length(union(genes_a,genes_b))
      
      return(j_ix)
    })
    
    j_ix_a = enframe(j_ix_a, "StudyB","JaccardIx") %>% 
      unnest()
    
    return(j_ix_a)
    
  }), "StudyA") %>% unnest()
  
  return(jaccard_res_sign)
  
}


























