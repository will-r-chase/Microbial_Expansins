library(purrr)
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(plyr)
library(dplyr)
library(magrittr)

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")
nodes_tree <- ggtree(tree) + geom_tiplab(size = 1) + geom_text2(aes(subset = !isTip, label = node), size = 1, hjust = -0.3)

#group HGT species to color branches
group_tree <- groupOTU(tree, .node=c("Ralstonia-syzygii", "Ralstonia-solanacearum", "Erwinia-tracheiphila", "Pantoea-stewartii", "Lonsdalea-quercina", "Pectobacterium-atrosepticum", "Pectobacterium-wasabiae", "Pectobacterium-parmentieri", "Pectobacterium-betavasculorum", "Pectobacterium-carotovorum", "Dickeya-zeae", "Dickeya-dadantii", "Dickeya-dianthicola", "Dickeya-chrysanthemi", "Dickeya-solani", "Cedecea-neteri", "Sphaerotilus-natas", "Janthinobacterium-sp", "Methylibium-sp", "Acidovorax-radicis", "Leptothrix-cholodnii", "Polyangium-brachysporum", "Vitrella-brassicaformis", "Sorangium-cellulosum", "Neocallimastix-californiae", "Anaeromyces-robustus", "Piromyces-finnis", "Allomyces-macrogynus", "Hamadaea-tsunoensis", "Streptomyces-acidiscabies", "Uliginosibacterium-gangwonense", "Haloferula-sp", "Physarum-polycephalum", "Acanthamoeba-castellanii", "Trichoderma-harzianum", "Trichoderma-virens", "Trichoderma-pseudokoningii", "Trichoderma-reesei", "Trichoderma-asperellum", "Trichoderma-atroviride", "Talaromyces-stipitatus2", "Penicillium-decumbens", "Penicillium-oxalicum", "Penicillium-brasilianum2", "Talaromyces-cellulolyticus3", "Talaromyces-marneffei", "Talaromyces-verruculosus2", "Thalassiosira-oceanica"))

#select nodes to collapse, label collapsed nodes with predominant species
nodes_to_collapse <- c(923, 938, 943, 951, 983, 990, 1001, 1019, 1115, 1134, 615, 788, 836, 846, 774)
labels = c("Xanthomonas spp", "Xanthomonas spp", "Paenibacillus spp", "Bacillus spp", "Bacillus spp", "Myxobacteria", "Myxobacteria", "Actinomycetes", "β-proteobacteria", "Actinomycetes", "Ascomycetes", "Amoebozoa", "Viridiplantae", "Stramenopiles", "Basidiomycetes")

#get number of children in each clade
tree_dat <- group_tree %>% 
  as.treedata() %>% 
  tidytree::as_data_frame()

children<-list()
for(i in 1:length(nodes_to_collapse)){
  children[[i]]<-offspring(tree_dat, .node = nodes_to_collapse[i])
}
names(children)<-labels
children <- map(children, ~.[-(which(grepl("/", .$label))), ])

number_collapsed <- map_int(children, ~nrow(.))

tree_labels <- paste(labels, paste0("(", paste(number_collapsed, "taxa)", sep = " ")), sep = " ")

#plot tree
collapse_tree<-
  ggtree(group_tree, aes(color=group)) + 
  geom_cladelabel(node = 903, label = "Prokaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 10, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = c("darkgrey", "red")) 

#collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- collapse(collapse_tree, node = nodes_to_collapse[i], clade_name=tree_labels[i])
}

#plot collapsed tree w/ points
p<-collapse_tree + 
  geom_point2(aes(subset = (node %in% nodes_to_collapse)), size = 6, shape = 23, fill = "blue") + 
  geom_nodelab(aes(subset = (node %in% nodes_to_collapse)), color = "black", hjust = -0.05, size = 6)
p
