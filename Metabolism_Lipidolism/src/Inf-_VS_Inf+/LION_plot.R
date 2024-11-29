# Title: LION plot
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
library(readxl)
library(ggplot2)
library(export) 
library(ggrepel)
library(tidyverse)
library(ggprism)
## Read data  ------------------------------------------------------------------
df_enrich_file <- read_csv("Human_DGCA_LION.csv")
## Enrichment plot   -----------------------------------------------------------
df_enrich <- df_enrich_file %>%
  filter(`FDR q-value` < 0.05) %>%
  mutate(value = -log(`FDR q-value`,base = 10)) %>%
  arrange(value) %>%
  mutate(Discription = factor(Discription,levels = Discription))

p1 <- ggplot(df_enrich,aes(x=Discription,y=value)) +
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_x_discrete(expand = c(0,0.8))  + # expand distance of y axis and first bar
  labs(x='',y='-LOG(FDR q-value)') +
  theme_prism()
pdf("Lion_enrich.pdf",
    width = 9,height = 9)
p1 
dev.off()
