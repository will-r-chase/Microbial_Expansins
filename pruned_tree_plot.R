library(purrr)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(plyr)
library(dplyr)
library(magrittr)
library(readxl)

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("350seqs_Linsi_rooted.newick")
nodes_tree <- ggtree(tree) + geom_tiplab(size = 1) + geom_text2(aes(label = node), size = 1, hjust = -0.3)

#read data, clean up columns, attach to tree
data <- read_xlsx("fileS4_microbe-data.xlsx")
data <- data[, 1:5]
colnames(data) <- c("organism", "group", "plant_path", "plant_associate", "ecology")

p <- ggtree(tree)
p <- p %<+% data

#group tips by taxonomy
taxonomy_tree <- split(p$data$label, p$data$group)
taxonomy_tree <- groupOTU(tree, taxonomy_tree, group_name = "taxonomy")

#define colors for taxonomy coloring
colors1<-c("#a6cee3", "#1f78b4", "#b15928", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#cccc3d")


ggtree(taxonomy_tree, aes(color=taxonomy)) + 
  geom_tiplab(color = "black", size = 0.8) +
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -30, x = 0) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black", size = 1) +
  scale_color_manual(values = colors1) +
  theme(legend.position = "right")

#group HGT species
HGT_tree <- groupOTU(taxonomy_tree, 
                          .node=c("Ralstonia-syzygii", "Ralstonia-solanacearum", "Erwinia-tracheiphila", "Pantoea-stewartii", "Lonsdalea-quercina", "Pectobacterium-atrosepticum", "Pectobacterium-wasabiae", "Pectobacterium-parmentieri", "Pectobacterium-betavasculorum", "Pectobacterium-carotovorum", "Dickeya-zeae", "Dickeya-dadantii", "Dickeya-dianthicola", "Dickeya-chrysanthemi", "Dickeya-solani", "Cedecea-neteri", "Sphaerotilus-natas", "Janthinobacterium-sp", "Methylibium-sp", "Acidovorax-radicis", "Leptothrix-cholodnii", "Polyangium-brachysporum", "Sorangium-cellulosum", "Neocallimastix-californiae", "Anaeromyces-robustus", "Piromyces-finnis", "Allomyces-macrogynus", "Hamadaea-tsunoensis", "Streptomyces-acidiscabies", "Uliginosibacterium-gangwonense", "Haloferula-sp", "Acanthamoeba-castellanii", "Trichoderma-harzianum", "Trichoderma-virens", "Trichoderma-pseudokoningii", "Trichoderma-reesei", "Trichoderma-asperellum", "Trichoderma-atroviride", "Talaromyces-stipitatus2", "Penicillium-decumbens", "Penicillium-oxalicum", "Penicillium-brasilianum2", "Talaromyces-cellulolyticus3", "Talaromyces-marneffei", "Talaromyces-verruculosus2"),
                          group_name = "HGT"
)

#plot tree
HGT_tree<-
  ggtree(HGT_tree, aes(color=taxonomy, size=HGT)) + 
  geom_treescale(width = 0.5, linesize = 1, fontsize = 5, y = -5, x = 0) +
  geom_nodelab(aes(label = label), color = "black", size = 1, hjust = -0.01) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values=c(0.7, 1.3)) + 
  geom_tiplab(color = "black", size = 1)

