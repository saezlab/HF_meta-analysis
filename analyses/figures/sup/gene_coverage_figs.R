# MIT License

# Copyright (c) [2020] [Jan D. Lanzer]
# jan.lanzer@biquant.uni-heidelberg.de

# Description: Plotting of gene coverage per study

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#loading prerequisites (output from samplesize_calculator.R and gene_coverage.R)
sample.size = readRDS(file ="data/clinical_description/sample_sizes.rds")
jaccard_df = readRDS(file = "data/figure_objects/jaccard_df_allgenes.rds") %>%
  mutate(value = ifelse(Var1 == Var2, 0, value))


##################### J'accard plot of gene coverage

plot.jaccard = ggplot(jaccard_df, aes(x = Var1, 
                                      y = Var2,
                                      fill = value)) +
  geom_tile() +
  theme_minimal() + xlab("Predicted") +
  ylab("Predictor") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        text = element_text(size =15))+
    scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, 'RdPu'),
                       limits=c(0, 1)) +
  labs(fill = "Jaccard \nIndex")

##################### Absolute Gene coverage per study

sample.size = sample.size %>% 
  mutate(study= factor(study,  levels = rev(as.character(sample.size$study))),
         genes= genes/1000)

# Data plotting
plot.coverage = ggplot(data = sample.size, 
                       aes(x= study, 
                           y= genes))+
  geom_col(fill = "#444C5C")+
  ylab(expression(Gene~coverage~(10^3))) + 
  theme_minimal()+
  scale_y_continuous(position = "right") + 
  coord_flip()+
  theme(axis.title.y = element_blank(),
        text = element_text(size =15),
        axis.title.x= element_text(size= 12)
  )


##################### Combine Coverage and Jaccard

plot.coverage.jaccard =  plot_grid(plot.coverage,plot.jaccard,
                                   ncol = 2,
                                   align = "h",
                                   labels = c("A","B"), 
                                   rel_widths = c(1,1.2))

pdf("data/figures/sup/SupplementalFigure2.pdf",
    width = 10,
    height= 5)
plot.coverage.jaccard
dev.off()

