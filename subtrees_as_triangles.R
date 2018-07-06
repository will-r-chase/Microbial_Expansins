library(purrr)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidytree)
library(magrittr)


##subtrees as triangles functions##
###########################################################################
get_offsprings <- function(node_to_collapse, phylo) {
  # todo: assert that node is scalar
  phylo %>%
    tidytree::as_data_frame() %>%
    tidytree::offspring(.node = node_to_collapse) %>%
    dplyr::pull(node)
}

get_offspring_tips <- function(phylo, node_to_collapse) {
  # todo: assert that node is scalar
  phylo %>%
    ggtree::fortify() %>%
    tidytree::offspring(.node = node_to_collapse) %>%
    dplyr::filter(isTip) %>%
    dplyr::pull(node)
}

remove_collapsed_nodes <- function(phylo, nodes_to_collapse) {
  nodes <- purrr::map(nodes_to_collapse, get_offsprings, phylo = phylo) %>% unlist()
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(!node %in% nodes)
}

get_collapsed_offspring_nodes_coordinates <- function(phylo, nodes) {
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(node %in% nodes) %>%
    dplyr::summarise(xmax = max(x), xmin = min(x), ymax = max(y), ymin = min(y))
}

get_collapsed_node_coordinates <- function(phylo, node_to_collapse) {
  # todo: assert that node is scalar
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(node == node_to_collapse) %>%
    dplyr::select(x, y)
}


get_triangle_coordinates_ <- function(node, phylo, mode = c("max", "min", "mixed")) {
  mode <- match.arg(mode)
  # for one
  tips_to_collapse <- get_offspring_tips(phylo, node)
  node_to_collapse_xy <- get_collapsed_node_coordinates(phylo, node)
  tips_to_collapse_xy <- get_collapsed_offspring_nodes_coordinates(phylo, tips_to_collapse)
  
  triange_df <- mode %>% switch(
    max = dplyr::data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmax, tips_to_collapse_xy$xmax),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    ),
    min = data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmin, tips_to_collapse_xy$xmin),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    ),
    mixed = data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmin, tips_to_collapse_xy$xmax),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    )
  )
  return(triange_df)
}

get_triangle_coordinates <- function(phylo, nodes, mode = c("max", "min", "mixed")) {
  mode <- match.arg(mode)
  # todo: make sure there is no conflict between nodes (nesting...)
  purrr::map(nodes, get_triangle_coordinates_, phylo = phylo, mode = mode) %>%
    dplyr::bind_rows(.id = "node_collapsed")
}
####################################################################################


##plotting tree##

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")
nodes_tree <- ggtree(tree) + geom_tiplab(size = 1) + geom_text2(aes(subset = !isTip, label = node), size = 1, hjust = -0.3)

#group HGT species to color branches
group_test<-groupOTU(tree, .node=c("Ralstonia-syzygii", "Ralstonia-solanacearum", "Erwinia-tracheiphila", "Pantoea-stewartii", "Lonsdalea-quercina", "Pectobacterium-atrosepticum", "Pectobacterium-wasabiae", "Pectobacterium-parmentieri", "Pectobacterium-betavasculorum", "Pectobacterium-carotovorum", "Dickeya-zeae", "Dickeya-dadantii", "Dickeya-dianthicola", "Dickeya-chrysanthemi", "Dickeya-solani", "Cedecea-neteri", "Sphaerotilus-natas", "Janthinobacterium-sp", "Methylibium-sp", "Acidovorax-radicis", "Leptothrix-cholodnii", "Polyangium-brachysporum", "Vitrella-brassicaformis", "Sorangium-cellulosum", "Neocallimastix-californiae", "Anaeromyces-robustus", "Piromyces-finnis", "Allomyces-macrogynus", "Hamadaea-tsunoensis", "Streptomyces-acidiscabies", "Uliginosibacterium-gangwonense", "Haloferula-sp", "Physarum-polycephalum", "Acanthamoeba-castellanii", "Trichoderma-harzianum", "Trichoderma-virens", "Trichoderma-pseudokoningii", "Trichoderma-reesei", "Trichoderma-asperellum", "Trichoderma-atroviride", "Talaromyces-stipitatus2", "Penicillium-decumbens", "Penicillium-oxalicum", "Penicillium-brasilianum2", "Talaromyces-cellulolyticus3", "Talaromyces-marneffei", "Talaromyces-verruculosus2", "Thalassiosira-oceanica"))

#select nodes to collapse and remove nodes from tree data
nodes_to_collapse <- c(923, 938, 943, 951, 981, 983, 990, 1001, 1019, 1115, 1134, 615, 788, 836, 846)
collapsed_tree_df <- group_test %>% 
  remove_collapsed_nodes(nodes = nodes_to_collapse)

#get coords for triangles
triangles_df <- group_test %>%
  get_triangle_coordinates(nodes_to_collapse, mode = "min")

#plot ggtree with collapsed nodes as triangles
triangles_plot <- 
  ggtree(collapsed_tree_df, aes(color=group)) +
  geom_cladelabel(node = 903, label = "Prokaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  geom_cladelabel(node = 610, label = "Eukaryotes", color = "black", barsize = 2, align = T, fontsize = 5, offset = 1) +
  xlim(0, 8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  geom_polygon(
    data = triangles_df,
    mapping = aes(group = node_collapsed),
    color = "#333333",
    fill = "steelblue"
  ) +
  theme(
    strip.background = element_blank()
  )



