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

data <- read_xlsx("D:/R projects/Microbial_Expansins/fileS4_microbe-data.xlsx")
data <- data %>%
  filter(Organism != "Vitrella-brassicaformis") %>%
  filter(Organism != "Emiliania-huxleyi") %>%
  select(1:5) 

colnames(data) <- c("organism", "group", "plant_path", "plant_associate", "ecology")



