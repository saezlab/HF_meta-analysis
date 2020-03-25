# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Description: Generates dictionary of studies/IDs
#' 
#' 

GEO_ID = c("GSE57345","GSE5406","Translatome","PRJNA477855",
           "PRJNA198165","PRJNA522417","GSE1869","GSE55296",
           "GSE42955","GSE26887","GSE3585","PRJNA291619",
           "PRJNA246308","GSE76701","GSE16499","GSE123976")

newID = c("Liu15_M","Hannenhalli06","vanHeesch19","Sweet18","Yang14",
          "Spurrell19","Kittleson05","Tarazon14","Molina-Navarro13","Greco12",     
          "Barth06","Schiano17","Liu15_R","Kim16","Kong10","Pepin19")

dictionary = data.frame(GEO_ID, newID,
                        stringsAsFactors = F)

save(dictionary, file = "data/dictionaryIDs.ro")
