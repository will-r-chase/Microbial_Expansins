library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(readxl)
library(purrr)

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

#read data, clean up columns, attach to tree
data <- read_xlsx("SupTab3_ecologyData.xlsx")

data_filtered <- data %>%
  select(1:6) %>%
  filter(Group != "Archaeplastids") %>%
  mutate(kingdom = case_when(Group2 == "Bacteria" ~ "bacteria",
                             Group2 == "Fungi" ~ "fungi",
                             Group2 == "Other" | Group2 == "Haptophyta" | Group2 == "Eukaryota" ~ "protists")) %>%
  filter(is.na(kingdom) == FALSE)

colnames(data_filtered) <- c("organism", "group", "group2", "plant_path", "plant_associated", "ecology", "kingdom")
data_filtered$plant_associated[which(data_filtered$plant_associated=="No" & data_filtered$plant_path=="Yes")] <- "Yes"
data_filtered$ecology[which(data_filtered$ecology=="Extremophile; Commensal")] <- "Extremophile"
data_filtered$ecology[data_filtered$organism=="Mastigocoleus-testarum"] <- "Marine"

#define colors for taxonomy coloring
colors <- c("#016c59", "#33a02c", "#091b77", "#cece3e", "#8c510a", "#c51b7d", "#252525", "#ff7f00", "#4b93c3", "#8861d4")

#make ecology counts table with colors defined
color_df <- data.frame(ecology = unique(data_filtered$ecology), color = "NA")
color_df <- color_df[order(color_df$ecology), ]
color_df$color <- colors

ecology_counts <- data_filtered %>% 
  truly_group_by(kingdom, ecology) %>%
  tally() %>%
  group_by(kingdom) %>%
  mutate(percent = round((n/sum(n))*100, digits = 1)) %>% 
  inner_join(., color_df, by="ecology")

ecology_split <- split(ecology_counts, ecology_counts$kingdom)

#make ecology breakdown barplots
ecology_plots <- lapply(ecology_split, function(x){
  ggplot(x, aes(x = ecology, y = n)) + 
    geom_col(fill = x$color) + 
    ylim(0, (max(x$n)*1.1)) +
    coord_flip() + 
    theme_light() +
    theme(legend.position = "right", panel.grid.major = element_blank()) +
    labs(x = "Ecological Niche", y = "Number of Species") +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = paste0(percent, "%"), hjust = - 0.05))
}
)  

ecology_plots[[1]]
ggsave("bacteria_ecology.tiff", dpi = 300)
ecology_plots[[2]]
ggsave("fungi_ecology.tiff", dpi = 300)
ecology_plots[[3]]
ggsave("protist_ecology.tiff", dpi = 300)

#make plant association barplots 
plant_associated_counts <- data_filtered %>% 
  filter(plant_associated != "UnkNown") %>%
  truly_group_by(kingdom, plant_associated) %>%
  tally() %>%
  group_by(kingdom) %>% 
  mutate(percent = round((n/sum(n))*100, digits = 1)) 

plant_associated_split <- split(plant_associated_counts, plant_associated_counts$kingdom)

plant_associated_plots <- lapply(plant_associated_split, function(x){
  ggplot(x, aes(x = kingdom, y = percent, fill = plant_associated, label = paste0(percent, "%"))) +
    scale_fill_manual(values = c("white", "springgreen4")) +
    geom_col(color = "black", size = 0.5) + 
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), legend.position = "none") +
    geom_text(position = position_stack(vjust = 0.5), size = 6)
}
)

plant_associated_plots[[1]]
ggsave("bacteria_associated.tiff", dpi = 300, units = "mm", width = 70, height = 170)
plant_associated_plots[[2]]
ggsave("fungi_associated.tiff", dpi = 300, units = "mm", width = 70, height = 170)
plant_associated_plots[[3]]
ggsave("protist_associated.tiff", dpi = 300, units = "mm", width = 70, height = 170)

