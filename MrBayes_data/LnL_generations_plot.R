library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(scales)

setwd("D:/R projects/Microbial Expansins/Microbial_Expansins/MrBayes_data")

run1 <- fread("10milGen_run1_noprune.p")
run2 <- fread("10milGen_run2_noprune.p")

combined <- data.table(gen = run1$Gen, run1 = run1$LnL, run2 = run2$LnL)

combined %>%
  filter(run1 > -110000) %>%
  ggplot(., aes(x = gen)) + 
  geom_line(aes(y = run1, color = "steelblue"), alpha = 0.5) + 
  geom_line(aes(y = run2, color = "orangered1"), alpha = 0.5) + 
  labs(x = "Generation", y = "Log Likelihood") +
  scale_color_discrete(name = "Runs", labels = c("Run 2", "Run 1")) +
  theme_light() +
  scale_x_continuous(labels = comma) +
  theme(legend.position = "right", axis.title = element_text(size = 16), legend.text = element_text(size = 12), legend.title = element_blank(), legend.key.size = unit(3, "line"))
