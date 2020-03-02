#  Copyright (c) [2019] [Ricardo O. Ramirez Flores, Jan Lanzer]
#  roramirezf@uni-heidelberg.de
#' In this script 14 independent data sets are combined in a single
#' object that will be used for the whole analysis in the paper
#' 
#' Disclaimer: Each dataset was processed in independent scripts, which
#' apply normalizations and annotations accordingly to the needs
#' 

source("src/data_utils.R") #general functions
source("src/misc_utils.R")

# 1. Creation of the list of data sets

METAheart = list() #Object initialization

##GSE76701
load(file = "data_processing/processed/GSE76701_counts.ro")
load(file = "data_processing/processed/GSE76701_targets.ro")

METAheart[["GSE76701"]] = list("GEX"= GSE76701_counts,
                                       "TARGETS" = GSE76701_targets)

#GSE57345: Gender/Age/Disease/HTx/DCM-ICM/from raw
load("data_processing/processed/GSE57345_targets.ro")
load("data_processing/processed/GSE57345_counts.ro")

METAheart[["GSE57345"]] = list("GEX" = GSE57345_counts, 
                               "TARGETS" = GSE57345_targets)

#GSE42955: Gender/DCM-ICM/HTx/from raw
load("data_processing/processed/GSE42955_counts.ro")
load("data_processing/processed/GSE42955_targets.ro")

METAheart[["GSE42955"]] = list("GEX" = GSE42955_counts, 
                               "TARGETS" = GSE42955_targets)

#GSE1869: DCM-ICM/HTx/LVAD from proccessed
#Jan version, less samples, from processed (not used)
#incomplete raw data
load("data_processing/processed/GSE1869_counts.ro")
load("data_processing/processed/GSE1869_targets.ro")

METAheart[["GSE1869"]] = list("GEX" = GSE1869_counts, 
                              "TARGETS" = GSE1869_targets)

#GSE3585: Jan manual targets file from paper/ data from raw *
# Gender/Age/
load("data_processing/processed/GSE3585_counts.ro")
load("data_processing/processed/GSE3585_target.ro")

METAheart[["GSE3585"]] = list("GEX" = GSE3585_counts, 
                              "TARGETS" = GSE3585_target)

METAheart[["GSE3585"]]$TARGETS = mutate(METAheart[["GSE3585"]]$TARGETS,
                                        "HTx" = "yes",
                                        "DCM" = ifelse(grepl("DCM",Diagnosis),
                                                      "yes","no"))

#GSE26887: Diabetes / Age (*Confounding) / Gender / from raw
load("data_processing/processed/GSE26887_targets.ro")
load("data_processing/processed/GSE26887_counts.ro")

METAheart[["GSE26887"]] = list("GEX" = GSE26887_counts, 
                               "TARGETS" = GSE26887_targets)

#GSE5406: 
load("data_processing/processed/GSE5406_targets.ro")
load("data_processing/processed/GSE5406_counts.ro")

METAheart[["GSE5406"]] = list("GEX" = GSE5406_counts, 
                              "TARGETS" = GSE5406_targets)

#GSE55296
load("data_processing/processed/GSE55296_targets.ro")
load("data_processing/processed/GSE55296_counts.ro")
colnames(GSE55296_targets)[1] = "Sample"

METAheart[["GSE55296"]] = list("GEX" = GSE55296_counts, 
                               "TARGETS" = GSE55296_targets)

METAheart[["GSE55296"]]$TARGETS = mutate(METAheart[["GSE55296"]]$TARGETS,
                                        "HTx" = "yes",
                                        "DCM" = ifelse(grepl("dilated",disease),
                                                       "yes","no"))

#PRNJA198165: Age / Gender, LVAD: All samples HF
load("data_processing/processed/PRJNA198165_target.ro")
load("data_processing/processed/PRJNA198165_count.ro")

METAheart[["PRJNA198165"]] = list("GEX" = PRJNA198165_count, 
                                  "TARGETS" = PRJNA198165_target)

#PRJNA246308 * : No age, very small data set, hard coded that column
load("data_processing/processed/PRJNA246308_target.ro")
load("data_processing/processed/PRJNA246308_count.ro")

METAheart[["PRJNA246308"]] = list("GEX" = PRJNA246308_count, 
                                  "TARGETS" = PRJNA246308_target[,-5])

#PRJNA291619 * :
load("data_processing/processed/PRJNA291619_target.ro")
load("data_processing/processed/PRJNA291619_count.ro")

METAheart[["PRJNA291619"]] = list("GEX" = PRJNA291619_count, 
                                  "TARGETS" = PRJNA291619_target)

#PRJNA477855 * Check name
load("data_processing/processed/PRJNA477855_target.ro")
load("data_processing/processed/PRJNA477855_count.ro")
METAheart[["PRJNA477855"]] = list("GEX" = PRJNA47855_count, 
                                  "TARGETS" = PRJNA477855_target)

#GSE126573/PRJNA522417: Explanted and with DCM
load("data_processing/processed/PRJNA522417_target_adult.ro")
load("data_processing/processed/PRJNA522417_count_adult.ro")

METAheart[["PRJNA522417"]] = list("GEX" = PRJNA522417_count_adult, 
                                  "TARGETS" = PRJNA522417_target_adult)

#Translatome *The same problem
load("data_processing/processed/Translatome_target.ro")
load("data_processing/processed/Translatome_count.ro")

METAheart[["Translatome"]] = list("GEX" = Translatome_count, 
                                  "TARGETS" = Translatome_targets)


# 2. Remove incomplete information : Only applies to translatome 

METAheart = lapply(METAheart, function(x){
  x$TARGETS = na.omit(x$TARGETS)
  return(x)
})

# 3. Order samples in GEX matrix based on target file

METAheart = lapply(METAheart, function(x){
  x$GEX = x$GEX[,x$TARGETS$Sample]
  return(x)
})

# 4. Generate QC boxplots

experiments = names(METAheart)

pdf(file = "other_res/boxplotsMETA.pdf",width = 13, height = 8)
for(x in experiments){
  boxplot(METAheart[[x]]$GEX, main=x)
}
dev.off()

# 5. Generate QC MDS

pdf(file = "other_res/MDSMETA.pdf",width = 6, height = 4)
for(x in experiments){
  myMDS(METAheart[[x]]$GEX,METAheart[[x]]$TARGETS, main=x)
}
dev.off()

# 6. Change IDs for the paper

load("data/dictionaryIDs.ro")

new_ids = as.data.frame(dictionary)
rownames(new_ids) = new_ids$GEO_ID

current_ids = names(METAheart)

names(METAheart) = new_ids[current_ids, "newID"]

# 7. Sample calculator

ct= c()
dcm= c()
icm = c()

for (study in names(METAheart)){
  ct= c(ct,dim(METAheart[[study]]$TARGETS %>% dplyr::filter(HeartFailure == "no"))[1])
  dcm = c(dcm,dim(METAheart[[study]]$TARGETS %>% dplyr::filter(DCM == "yes"))[1])
  icm = c(icm,dim(METAheart[[study]]$TARGETS %>% dplyr::filter(DCM == "no") %>% 
                    dplyr::filter(HeartFailure == "yes"))[1])
  
  ct1 = dim(METAheart[[study]]$TARGETS %>% filter(HeartFailure == "no"))[1]
}

sample.sizes = data.frame("study"= as.character(names(METAheart)), "CT"= ct, "DCM"= dcm, "ICM"= icm) %>% 
  mutate(total = rowSums(.[,2:4])) %>%
  arrange(desc(total))

print("total sample summary")
print(sample.sizes)

print("total samples")
print(sum(sample.sizes$total))

print("number of controls")
print(sum(sample.sizes$CT))

print("number of DCM")
print(sum(sample.sizes$DCM))

print("number of ICM")
print(sum(sample.sizes$ICM))


# 6. Save object 

saveRDS(METAheart, file = "data/METAheart.rds")



