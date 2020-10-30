library(tidyverse)

gs= read.csv(file = "data/meta_analysis_results_alamandadi.csv")



gs = as_tibble(gs) %>% rename(gene= gene_symbol) %>% select(gene)

write.csv(gs, "data/meta_alimdadai_geneset.csv", row.names = F)
write.csv(gs[1:500,], "data/meta_alimdadai_geneset500.csv", row.names = F)
