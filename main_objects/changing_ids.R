# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we map to definitive IDs for figures
#' 
#' Inputs :  METAheart object with all experiments unified in a list with target and GEX info
#' 
#' Outputs : METAheart mapped to def IDs
#' 

load("HGEX_data/dictionaryIDs.ro")

new_ids = as.data.frame(dictionary)
rownames(new_ids) = new_ids$GEO_ID

METAheart = readRDS(file = "HGEX_data/METAheart.rds") #main object

current_ids = names(METAheart)

names(METAheart) = new_ids[current_ids, "newID"]

saveRDS(METAheart, file = "HGEX_data/METAheart.rds")














 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
















