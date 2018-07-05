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

nodes_to_collapse <- c(923, 938)
collapsed_tree_df <- tree %>% 
  remove_collapsed_nodes(nodes = nodes_to_collapse)

triangles_df <- tree %>%
  get_triangle_coordinates(nodes_to_collapse, mode = "mixed")

ggtree(collapsed_tree_df) +
  geom_polygon(
    data = triangles_df,
    mapping = aes(group = node_collapsed, fill = node_collapsed),
    color = "#333333"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    strip.background = element_blank()
  )
