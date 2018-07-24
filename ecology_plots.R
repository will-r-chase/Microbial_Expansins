library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ape)
library(ggtree)
library(readxl)
library(purrr)

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
p <- ggtree(tree)
p <- p %<+% data

#define colors for taxonomy coloring
colors <- c("#016c59", "#33a02c", "#e31a1c", "#6a3d9a", "#FFFF4D", "#8c510a", "#c51b7d", "#252525", "#ff7f00", "#1f78b4")

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
    coord_flip() + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank()) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(percent, "%"), hjust = - 0.05))
  }
)  

#group tips by taxonomy
eco_tree <- split(p$data$label, p$data$ecology)
eco_tree <- groupOTU(tree, eco_tree, group_name = "ecology")

#plot tree
ecology_tree<-
  ggtree(eco_tree, aes(color=ecology)) + 
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>=70 & as.numeric(sub(".*/", "", label))>=95 & !isTip), color = "black") +
  xlim(0, 8) +
  scale_color_manual(values = colors) +
  theme(legend.position = "right") +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) +
  geom_strip(1, 1, label = "", barsize = 2, color = "black", align = T, fontsize = 10, offset = 1) 

