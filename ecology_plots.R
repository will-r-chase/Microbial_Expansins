library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ape)
library(ggtree)
library(readxl)
library(purrr)
library(gridExtra)

##function to fix zero-length groups in ecology plots
truly_group_by <- function(data, ...){
  dots <- quos(...)
  data <- group_by( data, !!!dots )
  
  labels <- attr( data, "labels" )
  labnames <- names(labels)
  labels <- mutate( labels, ..index.. =  attr(data, "indices") )
  
  expanded <- labels %>%
    tidyr::expand( !!!dots ) %>%
    left_join( labels, by = labnames ) %>%
    mutate( ..index.. = map(..index.., ~if(is.null(.x)) integer() else .x ) )
  
  indices <- pull( expanded, ..index..)
  group_sizes <- map_int( indices, length)
  labels <- select( expanded, -..index..)
  
  attr(data, "labels")  <- labels
  attr(data, "indices") <- indices
  attr(data, "group_sizes") <- group_sizes
  
  data
}

#read tree and find nodes to collapse, only monophyletic clades can be collapsed
tree <- read.tree("fileS1_msa_fixed_rooted.tree")

#read data, clean up columns, attach to tree
data <- read_xlsx("fileS4_microbe-data.xlsx")
data <- data[, 1:5]
colnames(data) <- c("organism", "group", "plant_path", "plant_associate", "ecology")
data$plant_associate[which(data$plant_associate=="No" & data$plant_path=="Yes")] <- "Yes"
p <- ggtree(tree)
p <- p %<+% data

#define colors for taxonomy coloring
colors <- c("#016c59", "#33a02c", "#e31a1c", "#6a3d9a", "#cece3e", "#8c510a", "#c51b7d", "#252525", "#ff7f00", "#1f78b4")

#make ecology counts table with colors defined
color_df <- data.frame(ecology = unique(data$ecology), color = "NA")
color_df <- color_df[order(color_df$ecology), ]
color_df$color <- colors

ecology_counts <- data %>% 
  truly_group_by(group, ecology) %>%
  tally() %>%
  group_by(group) %>% 
  mutate(percent = round((n/sum(n))*100, digits = 1)) %>% 
  inner_join(., color_df, by="ecology")

ecology_split <- split(ecology_counts, ecology_counts$group)

ecology_plots <- lapply(ecology_split, function(x){
  ggplot(x, aes(x = ecology, y = n)) + 
    geom_bar(stat = "identity", fill=x$color) + 
    expand_limits(y = max(x$n)*1.2) +
    coord_flip() + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(percent, "%"), hjust = - 0.05))
  }
)  

#group tips by taxonomy
eco_split <- split(p$data$label, p$data$ecology)
plant_associate_split <- split(p$data$label, p$data$plant_associate)
plant_path_split <- split(p$data$label, p$data$plant_path)
eco_tree <- groupOTU(tree, eco_split, group_name = "ecology")
eco_tree <- groupOTU(eco_tree, plant_associate_split, group_name = "plant_associate")
eco_tree <- groupOTU(eco_tree, plant_path_split, group_name = "plant_path")

#plot tree w/ legend
ecology_tree <-
  ggtree(eco_tree, aes(color=ecology)) + 
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black", size = 0.8) +
  geom_tippoint(aes(shape = plant_associate), size = 0.8, color = "seagreen3") +
  scale_shape_manual(values = c(32, 17)) +
  xlim(0, 8) +
  scale_color_manual(values = colors) +
  theme(legend.position = "right") +
  geom_strip(296, 315, label = "Xanthomonads", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(317, 347, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(349, 360, label = "Enterobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(363, 374, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(380, 394, label = "Myxobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(426, 500, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(527, 602, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(2, 161, label = "Ascomycetes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(163, 174, label = "Basidiomycetes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(178, 207, label = "Amoebozoa", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(225, 231, label = "Archaeplastida", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(259, 291, label = "Stramenopiles", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1)

#plot tree w/o legend
ecology_tree2 <-
  ggtree(eco_tree, aes(color=ecology)) + 
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black") +
  geom_tippoint(aes(shape = plant_associate), size = 0.8, color = "seagreen3") +
  scale_shape_manual(values = c(32, 17)) +
  xlim(0, 8) +
  scale_color_manual(values = colors) +
  theme_tree2() +
  geom_strip(296, 315, label = "Xanthomonads", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(317, 347, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(349, 360, label = "Enterobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(363, 374, label = "Firmicutes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(380, 394, label = "Myxobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(426, 500, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(527, 602, label = "Actinobacteria", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(2, 161, label = "Ascomycetes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(163, 174, label = "Basidiomycetes", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(178, 207, label = "Amoebozoa", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(225, 231, label = "Archaeplastida", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(259, 291, label = "Stramenopiles", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1)

#setup grid layout for multiplot
lay <- rbind(
             c(1,1,1,1,2),
             c(1,1,1,1,3),
             c(1,1,1,1,4),
             c(1,1,1,1,5),
             c(1,1,1,1,6),
             c(1,1,1,1,7),
             c(1,1,1,1,8),
             c(1,1,1,1,9),
             c(1,1,1,1,10),
             c(1,1,1,1,11),
             c(1,1,1,1,12),
             c(1,1,1,1,13)
             )

#reorder grobs for multiplot
order <- c("tree", "Xanthomonads", "Firmicutes", "Enterobacteria", "Myxobacteria", "Actinobacteria", "Other", "Ascomycetes", "Basidiomycetes", "Amoebozoa", "Archaeplastids", "Stramenopiles")

#plot grobs
plots_list <- ecology_plots
plots_list$tree <- ecology_tree2
plots_list <- plots_list[order]
grid.arrange(grobs = plots_list, layout_matrix = lay)



##testing insets... not working##

#simple version
ecology_plots_simple <- lapply(ecology_split, function(x){
  ggplot(x, aes(x = ecology, y = n)) + 
    geom_bar(stat = "identity") + 
    theme_inset()
}
)  

#get nodes to place bar chart insets
names(ecology_plots_simple) <- c(1025, 788, 837, 616, 1125, 774, 969, 916, 1001, 1117, 849, 921)
test<-ggtree(tree)
inset_tree <- inset(test, ecology_plots_simple, width = 100, height = 100, x = "branch")

#troubleshooting insets
tr <- rtree(15)
v <- ggtree(tr)
d <- lapply(1:15, rnorm, n=100)
ylim <- range(unlist(d))
bx <- lapply(d, function(y) {
  dd <- data.frame(y=y)
  ggplot(dd, aes(x=1, y=y))+geom_boxplot() + ylim(ylim) + theme_inset()
})
names(bx) <- 1:15
inset(test, bx, width=1, height=1)
