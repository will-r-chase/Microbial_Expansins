##plotting tree##

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")
nodes_tree <- ggtree(tree) + geom_tiplab(size = 1) + geom_text2(aes(subset = !isTip, label = node), size = 1, hjust = -0.3)

#group HGT species to color branches
group_tree <- groupOTU(tree, .node=c("Ralstonia-syzygii", "Ralstonia-solanacearum", "Erwinia-tracheiphila", "Pantoea-stewartii", "Lonsdalea-quercina", "Pectobacterium-atrosepticum", "Pectobacterium-wasabiae", "Pectobacterium-parmentieri", "Pectobacterium-betavasculorum", "Pectobacterium-carotovorum", "Dickeya-zeae", "Dickeya-dadantii", "Dickeya-dianthicola", "Dickeya-chrysanthemi", "Dickeya-solani", "Cedecea-neteri", "Sphaerotilus-natas", "Janthinobacterium-sp", "Methylibium-sp", "Acidovorax-radicis", "Leptothrix-cholodnii", "Polyangium-brachysporum", "Vitrella-brassicaformis", "Sorangium-cellulosum", "Neocallimastix-californiae", "Anaeromyces-robustus", "Piromyces-finnis", "Allomyces-macrogynus", "Hamadaea-tsunoensis", "Streptomyces-acidiscabies", "Uliginosibacterium-gangwonense", "Haloferula-sp", "Physarum-polycephalum", "Acanthamoeba-castellanii", "Trichoderma-harzianum", "Trichoderma-virens", "Trichoderma-pseudokoningii", "Trichoderma-reesei", "Trichoderma-asperellum", "Trichoderma-atroviride", "Talaromyces-stipitatus2", "Penicillium-decumbens", "Penicillium-oxalicum", "Penicillium-brasilianum2", "Talaromyces-cellulolyticus3", "Talaromyces-marneffei", "Talaromyces-verruculosus2", "Thalassiosira-oceanica"))

#select nodes to collapse 
nodes_to_collapse <- c(923, 938, 943, 951, 981, 983, 990, 1001, 1019, 1115, 1134, 615, 788, 836, 846, 774)

#plot tree
collapse_tree<-
  ggtree(group_test, aes(color=group)) + 
  geom_cladelabel(node = 903, label = "Prokaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonas spp", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 938), label = "Xanthomonas spp", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 943), label = "Paenibacillus spp", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +
  geom_nodelab(aes(subset = node == 923), label = "Xanthomonads", color = "black", hjust = -0.4) +

collapse_tree<-
  ggtree(group_test, aes(color=group)) + 
  geom_cladelabel(node = 903, label = "Prokaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  geom_tiplab(color = "black")

#collapse nodes
for(i in 1:length(nodes_to_collapse)){
  collapse_tree <- ggtree::collapse(collapse_tree, node = nodes_to_collapse[i])
}

#plot collapsed tree w/ points
p<-collapse_tree + geom_point2(aes(subset = (node %in% nodes_to_collapse)), size = 3, shape = 23, fill = "blue")

labels = LETTERS[1:16]

p + geom_nodelab2(aes(subset = (node %in% nodes_to_collapse), label = labels))


#get number of children in each clade
tree_dat <- group_test %>% 
  as.treedata() %>% 
  tidytree::as_data_frame()

children<-list()

for(i in 1:length(nodes_to_collapse)){
  children2[[i]]<-offspring(tree_dat, .node = nodes_to_collapse[i])
}

test_df<-children2[[2]]
test_df<-test_df[-(which(grepl("/", test_df$label))), ]

names(children2) <- c("Xanthomonads", "Xanthomonads 2", "Paenibacillus", )
