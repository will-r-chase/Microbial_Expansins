library(treeio)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(readxl)
library(dplyr)

setwd("D:/R projects/Microbial_Expansins/MrBayes_data")

bayes <- read.mrbayes("10milGen_noprune_contree_rooted.nexus")
bayes
nodes_tree <- ggtree(bayes) + geom_tiplab(size = 1) + geom_text2(aes(label = node), size = 1, hjust = -0.3)
nodes_tree

data <- read_xlsx("D:/R projects/Microbial Expansins/Microbial_Expansins/fileS4_microbe-data.xlsx")
data <- data %>%
  filter(Organism != "Vitrella-brassicaformis") %>%
  filter(Organism != "Emiliania-huxleyi") %>%
  select(1:5) 

colnames(data) <- c("organism", "group", "plant_path", "plant_associate", "ecology")

bayes_data<-bayes %>% tidytree::as_data_frame()

data2 <- data[which(data$organism%in%bayes_data$label),]		##remove any names in your data frame not found in tip labels (eg. node names)
tr_tips <- data.frame(organism = as.character(bayes_data$label))  ##make data frame with tree tip labels
data3 <- left_join(tr_tips, data2, by="organism")		##reorder your data frame to match order or tip labels

tree_split <- split(bayes_data$label, data3$group)
tax_group_tree <- groupOTU(bayes_data, tree_split, group_name = "taxonomy")

tax_group_tree$upper_prob <- tax_group_tree$prob_range %>%
  purrr::map(., 2, .null = NA) %>%
  unlist() %>%
  round(., digits = 2)

tax_treedata <- as.treedata(tax_group_tree)

colors1<-c("#a6cee3", "#1f78b4", "#b15928", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#cccc3d")
colors2<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#FFFF4D", "#b15928")

ggtree(tax_treedata, aes(color = taxonomy)) + 
  geom_tiplab(size = 1, color = "black") + 
  geom_nodelab(aes(label = prob_percent), color = "black", size = 1.5, hjust = 0.1) +
  geom_treescale(y = -30, width = 0.5) +
  theme(legend.position = "right") +
  scale_color_manual(values = colors1)
