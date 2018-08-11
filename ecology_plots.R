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
colnames(data) <- c("organism", "group", "plant_path", "plant_associated", "ecology")
data$plant_associated[which(data$plant_associated=="No" & data$plant_path=="Yes")] <- "Yes"
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

#make ecology breakdown barplots
ecology_plots <- lapply(ecology_split, function(x){
  ggplot(x, aes(x = ecology, y = percent)) + 
    geom_col(fill = x$color) + 
    ylim(0, 115) +
    coord_flip() + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(percent, "%"), hjust = - 0.05))
  }
)  

#make plant association barplots 
plant_associated_counts <- data %>% 
  truly_group_by(group, plant_associated) %>%
  tally() %>%
  group_by(group) %>% 
  mutate(percent = round((n/sum(n))*100, digits = 1)) 

plant_associated_split <- split(plant_associated_counts, plant_associated_counts$group)

plant_associated_plots <- lapply(plant_associated_split, function(x){
  ggplot(x, aes(x = group, y = percent, fill = plant_associated)) +
    scale_fill_manual(values = c("white", "springgreen4")) +
    geom_col(color = "black", size = 0.5) + 
    theme_inset() +
    labs(y = paste0(x$percent[2], "%"))
  }
)

plant_associated_plots <- plant_associated_plots[-c(5, 10)]
names(plant_associated_plots) <- c(426, 178, 225, 2, 163, 349, 317, 380, 259, 296)

#group tips by taxonomy
eco_split <- split(p$data$label, p$data$ecology)
plant_associated_split <- split(p$data$label, p$data$plant_associated)
plant_path_split <- split(p$data$label, p$data$plant_path)
eco_tree <- groupOTU(tree, eco_split, group_name = "ecology")
eco_tree <- groupOTU(eco_tree, plant_associated_split, group_name = "plant_associated")
eco_tree <- groupOTU(eco_tree, plant_path_split, group_name = "plant_path")

#plot tree w/ legend
ecology_tree <-
  ggtree(eco_tree, aes(color=ecology)) + 
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black", size = 0.8) +
  geom_tippoint(aes(shape = plant_associated), size = 0.8, color = "seagreen3") +
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
  xlim(0, 6) +
  scale_color_manual(values = colors) +
  geom_treescale(width = 0.5, linesize = 3, fontsize = 5, y = 0, x = 0) +
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

#plot plant_association as inset
ecology_tree_insets <- inset(ecology_tree2, plant_associated_plots, width = 0.2, height = 30, hjust = -1.2)

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
             c(1,1,1,1,12)
             )

#reorder grobs for multiplot
order <- c("tree", "Xanthomonads", "Firmicutes", "Enterobacteria", "Myxobacteria", "Actinobacteria", "Ascomycetes", "Basidiomycetes", "Amoebozoa", "Archaeplastids", "Stramenopiles")

#plot grobs
plots_list <- ecology_plots
plots_list$tree <- ecology_tree_insets
plots_list <- plots_list[order]
grid.arrange(grobs = plots_list, layout_matrix = lay)
