# Author: Jan Lanzer, 2019
# Description: Gathering of clinical information (Age, Gender, EF) f HF Studies

## load packages and objects ##
library(tidyverse)
library(ggpubr)
library(ggforce)
METAheart = readRDS(file = "data/METAheart.rds")
clinicalboolean = as.tibble(read.table("data/clinical_description/Tables for MetaHeart Manuscript - Table 2 - Clinical Characteristics.csv",
                                       header = T,
                                       sep = ",",
                                       stringsAsFactors = F)) 
load("data/dictionaryIDs.ro")
dictionary = dictionary  %>% as.data.frame() %>% column_to_rownames("newID")

names(METAheart) = dictionary[names(METAheart),] # rename metaehart object to GEO IDs(old version)
dictionary = dictionary  %>% rownames_to_column("newID") 
Experiments = names(METAheart)

######################## Transforming Clinicalboolean Matrix
#rename the studies based on the IDs from dictionary
clinicalboolean = clinicalboolean %>% 
  inner_join(.,dictionary %>% rename(Study.ID = GEO_ID)) %>% 
  select(-Study.ID) %>% 
  rename(Study.ID= newID) %>% 
  mutate(Study.ID = factor(Study.ID, levels = c(rev(dictionary$newID))))

######################## Creating Age Data Frame ######W##################

# Ages is a table with information of the means and sd of CT & HF samples for each study
ages = as.tibble(data.frame(Study= names(METAheart))) %>% 
  mutate(HFMeanAge= 0,
  HFSDAge= 0,
  NFMeanAge= 0,
  NFSDAge= 0) %>% 
  column_to_rownames("Study")


# adding age information for Liu_R that was excluded by Rico in the Metaheart object
# because it is not suitable for correction in limma
load("data/PRNJA246308_target.ro")

METAheart$PRJNA246308$TARGETS = METAheart$PRJNA246308$TARGETS %>% 
  left_join(PRNJA246308_target %>% select(Sample, Age), by= "Sample")


# For those studies, where samoplewise age information is available, the summary statistics will
# be calculated in the following for-loop.

for (name in names(METAheart)){
  ages[name,"HFMeanAge"]= mean(METAheart[[name]]$TARGETS$Age[METAheart[[name]]$TARGETS$HeartFailure == "yes"], na.rm =TRUE )
  ages[name,"HFSDAge"]= sd(METAheart[[name]]$TARGETS$Age[METAheart[[name]]$TARGETS$HeartFailure == "yes"], na.rm =TRUE )
  ages[name,"NFMeanAge"]= mean(METAheart[[name]]$TARGETS$Age[METAheart[[name]]$TARGETS$HeartFailure == "no"], na.rm =TRUE )
  ages[name,"NFSDAge"]= sd(METAheart[[name]]$TARGETS$Age[METAheart[[name]]$TARGETS$HeartFailure == "no"], na.rm =TRUE )
}


# Summary statistics that were manually curated from publications will be added. 
#GSE5406
ages["GSE5406",] = c(54,12,57,12)
#GSE55296
ages["GSE55296",] = c(52.5,9.16,47,16)
#PRJNA291619
ages["PRJNA291619",] = c(57.5, 10.2, 57.5, 11.7)
#GSE42955
ages["GSE42955",] = c(NA,NA,55,3)
#GSE16499
ages["GSE16499",] = c(52, 4, 49,8 )
#GSE1869
#ages["GSE1869","NFMeanAge"] = mean(GSE1869_target$Age[GSE1869_target$HeartFailure== "no"], na.rm = T)
#ages["GSE1869","NFSDAge"] = sd(GSE1869_target$Age[GSE1869_target$HeartFailure== "no"], na.rm = T)

#Round age
ages = round(ages, 0)

# Save ages as R Object
dictionary = dictionary %>% column_to_rownames("GEO_ID")
rownames(ages) = dictionary[rownames(ages), "newID"]

save(ages, file = "data/clinical_description/age_matrix.ro")

######################## Creating Gender Data Frame ########################
# gender is a table with information of the percentage of female samples for each study
gender = as.tibble(data.frame(Study= Experiments))%>% 
  mutate(female_HF= NA) %>%
  mutate(female_NF= NA) %>%
  column_to_rownames("Study")

# For those studies, where samoplewise gender information is available, the proportion of females will 
# be calculated in the following for-loop.
for (name in Experiments){
    # Caluclate gender of diseased groups
  hfgender = table( METAheart[[name]]$TARGETS$Gender[METAheart[[name]]$TARGETS$HeartFailure == "yes"])
    hffemales = hfgender["female"]*100 / (hfgender["female"]+ hfgender["male"])  
  gender[name,1] = hffemales
  
  # Calculate gender of control groups
  nfgender = table(METAheart[[name]]$TARGETS$Gender[METAheart[[name]]$TARGETS$HeartFailure == "no"])
  nffemales = nfgender["female"]*100 / (nfgender["female"]+ nfgender["male"])  
  gender[name,2] = nffemales 
  
}

# Round proportions
gender = round(gender, 0)

# Manually add information curated from papers (summary statistics)
gender["GSE55296", ] = c(4,20)
gender["PRJNA291619",] = c(50,75)
gender["GSE1869",] = c(34,NA)
gender["GSE42955",] = c(0,20) # added manually because all samples were males 
gender["GSE16499", ] =c(20,20)
gender["GSE123976",] =c(0,67)

# Save gender as R Object
rownames(gender) = dictionary[rownames(gender), "newID"]
save(gender, file = "data/clinical_description/gender_matrix.ro")
