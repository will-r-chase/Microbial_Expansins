library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)

colnames(data) <- c("organism", "group", "kingdom", "plant_pathogen", "plant_associated", "ecology") 
data <- data[, 1:6]

ecology_counts <- data %>% 
  group_by(group, ecology) %>% 
  tally() %>% 
  group_by(group) %>% 
  mutate(percent = round((n/sum(n))*100, digits = 1)) 
 
ecology_split <- split(ecology_counts, ecology_counts$group)

ecology_plots <- lapply(ecology_split, function(x){
  ggplot(x, aes(x = reorder(ecology, n), y = n, fill = ecology)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_few() + 
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Ecology", y = "Count") + 
    geom_text(aes(label = paste0(percent, "%"), hjust = - 0.05))
  }
)  

arranged_ecology_plots<-marrangeGrob(grobs=ecology_plots, nrow=3, ncol=3)
ggsave("ecology_plots.pdf", arranged_ecology_plots, height=20, width=20, units="in")


