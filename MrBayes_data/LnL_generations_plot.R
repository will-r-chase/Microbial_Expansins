library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)

setwd("D:/R projects/Microbial Expansins/Microbial_Expansins/MrBayes_data")

run1 <- fread("10milGen_run1_noprune.p")
run2 <- fread("10milGen_run2_noprune.p")

combined <- data.table(gen = run1$Gen, run1 = run1$LnL, run2 = run2$LnL)

combined %>%
  filter(run1 > -110000) %>%
  ggplot(., aes(x = gen)) + 
  geom_line(aes(y = run1, color = "steelblue"), alpha = 0.5) + 
  geom_line(aes(y = run2, color = "orangered1"), alpha = 0.5) + 
  labs(x = "Generation", y = "-LnL") +
  scale_color_discrete(name = "Runs", labels = c("Run 2", "Run 1")) +
  theme_bw() +
  theme(legend.position = "right")
