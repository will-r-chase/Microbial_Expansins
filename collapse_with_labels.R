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
tree <- read.tree("fileS1_msa_fixed_rooted.tree")
nodes_tree <- ggtree(tree) + geom_tiplab(size = 1) + geom_text2(aes(label = node), size = 1, hjust = -0.3)
nodes_to_collapse <- c(923, 938, 943, 951, 983, 990, 1001, 1019, 1115, 1134, 615, 788, 836, 846, 774)
labels = c("Xanthomonas spp", "Xanthomonas spp", "Paenibacillus spp", "Bacillus spp", "Bacillus spp", "Myxobacteria", "Myxobacteria", "Actinobacteria", "β-proteobacteria", "Actinobacteria", "Ascomycetes", "Amoebozoa", "Viridiplantae", "Stramenopiles", "Basidiomycetes")

#read data, clean up columns, attach to tree
data <- read_xlsx("fileS4_microbe-data.xlsx")
data <- data[, 1:5]
colnames(data) <- c("organism", "group", "plant_path", "plant_associate", "ecology")
p <- ggtree(tree)
p <- p %<+% data

#group tips by taxonomy
taxonomy_tree <- split(p$data$label, p$data$group)
taxonomy_tree <- groupOTU(tree, taxonomy_tree, group_name = "taxonomy")

#group HGT species
taxonomy_tree <- groupOTU(taxonomy_tree, 
                       .node=c("Ralstonia-syzygii", "Ralstonia-solanacearum", "Erwinia-tracheiphila", "Pantoea-stewartii", "Lonsdalea-quercina", "Pectobacterium-atrosepticum", "Pectobacterium-wasabiae", "Pectobacterium-parmentieri", "Pectobacterium-betavasculorum", "Pectobacterium-carotovorum", "Dickeya-zeae", "Dickeya-dadantii", "Dickeya-dianthicola", "Dickeya-chrysanthemi", "Dickeya-solani", "Cedecea-neteri", "Sphaerotilus-natas", "Janthinobacterium-sp", "Methylibium-sp", "Acidovorax-radicis", "Leptothrix-cholodnii", "Polyangium-brachysporum", "Vitrella-brassicaformis", "Sorangium-cellulosum", "Neocallimastix-californiae", "Anaeromyces-robustus", "Piromyces-finnis", "Allomyces-macrogynus", "Hamadaea-tsunoensis", "Streptomyces-acidiscabies", "Uliginosibacterium-gangwonense", "Haloferula-sp", "Physarum-polycephalum", "Acanthamoeba-castellanii", "Trichoderma-harzianum", "Trichoderma-virens", "Trichoderma-pseudokoningii", "Trichoderma-reesei", "Trichoderma-asperellum", "Trichoderma-atroviride", "Talaromyces-stipitatus2", "Penicillium-decumbens", "Penicillium-oxalicum", "Penicillium-brasilianum2", "Talaromyces-cellulolyticus3", "Talaromyces-marneffei", "Talaromyces-verruculosus2", "Thalassiosira-oceanica"),
                       group_name = "HGT"
                       )

#get number of children in each clade
tree_dat <- taxonomy_tree %>% 
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

#define colors for taxonomy coloring
colors1<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#cece3e", "#b15928")
colors2<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#FFFF4D", "#b15928")

#plot tree
collapse_tree<-
  ggtree(taxonomy_tree, aes(color=taxonomy, size=HGT)) + 
  geom_treescale(width = 0.5, linesize = 2, fontsize = 5, y = 0, x = 0) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black") +
  geom_cladelabel(node = 903, label = "Prokaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  geom_strip(604, 608, label = "Mixed Prokaryotes/Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values=c(1, 3)) + 
  geom_tiplab(color = "black", size = 3)

#collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- collapse(collapse_tree, node = nodes_to_collapse[i], clade_name=tree_labels[i])
}

#fix wrongly labeled clades
collapse_tree$data$taxonomy[990]<-"Myxobacteria"
collapse_tree$data$taxonomy[1019]<-"Actinobacteria"
collapse_tree$data$taxonomy[1115]<-"B-proteobacteria"

#plot collapsed tree w/ points
p2<-collapse_tree + 
  geom_point2(aes(subset = (node %in% nodes_to_collapse), fill=taxonomy), size = 6, shape = 23) + 
  geom_nodelab(aes(subset = (node %in% nodes_to_collapse)), color = "black", hjust = -0.05, size = 6) +
  scale_fill_manual(values = colors2)
p2


ggtree(taxonomy_tree, aes(color=taxonomy)) + geom_tiplab(size=1, hjust = -0.3) + scale_color_manual(values = colors1) + geom_text2(aes(label = node), size = 1, hjust = -0.3, color = "black") + theme(legend.position = "right")
