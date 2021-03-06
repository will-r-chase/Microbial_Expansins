library(tidyverse)
library(ape)
library(tidytree)
library(ggtree)
library(readxl)

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")
nodes_to_collapse <- c(943, 951, 968, 1001, 1135, 616, 788, 846, 1022)
labels = c("Paenibacillus", "Bacillus", "Enterobacteria", "Myxobacteria", "Actinobacteria", "Ascomycetes", "Amoebozoa", "Stramenopiles", "Actinobacteria")

#read data, clean up columns, attach to tree
data_fusion <- read_xlsx("fusions_for_tree.xlsx")
p <- ggtree(tree)
p <- p %<+% data_fusion

data_tax <- read_xlsx("fileS4_microbe-data.xlsx")
data_tax <- data_tax[, 1:5]
colnames(data_tax) <- c("organism", "group", "plant_path", "plant_associate", "ecology")
p <- p %<+% data_tax

#group tips by taxonomy
taxonomy_tree <- split(p$data$label, p$data$group)
taxonomy_tree <- groupOTU(tree, taxonomy_tree, group_name = "taxonomy")

#group tips by fusion
fusion_tree <- split(p$data$label, p$data$`Domain Architecture`)
fusion_tree <- groupOTU(taxonomy_tree, fusion_tree, group_name = "fusion")

#define colors for taxonomy coloring
colors1<-c("#a6cee3", "#1f78b4", "#b15928", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#cccc3d")
colors2<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#FFFF4D", "#b15928")
shapes <- c(32, 17, 17, 17, 17, 17, 17, 17)

treedat <- tidytree::as_data_frame(fusion_tree)

fusion_nodes <-
  treedat %>%
  filter(fusion != 0) %>%
  pull(node)

normal_nodes <-
  treedat %>%
  filter(fusion == 0) %>%
  pull(node)


###this is for collapsing, edit later for fusions
tree_dat <- fusion_tree %>% 
  as.treedata() %>% 
  tidytree::as_data_frame()

children<-list()

for(i in 1:length(nodes_to_collapse)){
  children[[i]]<-offspring(tree_dat, .node = nodes_to_collapse[i])
}

names(children)<-labels

children <- map(children, ~.[-(which(grepl("/", .$label))), ])

number_collapsed <- map_int(children, ~nrow(.))

#make labels for collapsed clades
tree_labels <- paste(labels, paste0("(", paste(number_collapsed, "taxa)", sep = " ")), sep = " ")


collapse_tree <- 
  ggtree(fusion_tree, aes(color=taxonomy)) + 
  geom_tiplab(aes(subset = node %in% fusion_nodes & isTip, label = paste(label, fusion, sep = ", ")), color = "red", size = 2, fontface = "bold", hjust = -0.05) +
  geom_tiplab(aes(subset = node %in% normal_nodes & isTip, label = label), color = "black", size = 2, hjust = -0.05) +
  #geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -30, x = 0) +
  geom_point2(aes(subset = isTip, shape = fusion), color = "red", size = 2) +
  scale_color_manual(values = colors1) +
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -5, x = 0) +
  scale_shape_manual(values = shapes) +
  xlim(0, 8) +
  theme(legend.position = "right")

#collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- collapse(collapse_tree, node = nodes_to_collapse[i], clade_name=tree_labels[i])
}

#fix wrongly labeled clades
# collapse_tree$data$taxonomy[990]<-"Myxobacteria"
# collapse_tree$data$taxonomy[1019]<-"Actinobacteria"
# collapse_tree$data$taxonomy[1115]<-"B-proteobacteria"

#plot collapsed tree w/ points
p2<-collapse_tree + 
  geom_point2(aes(subset = (node %in% nodes_to_collapse), fill=taxonomy), size = 3, shape = 23) + 
  geom_nodelab(aes(subset = (node %in% nodes_to_collapse)), color = "black", hjust = -0.05, size = 3) +
  scale_fill_manual(values = colors2)
p2
