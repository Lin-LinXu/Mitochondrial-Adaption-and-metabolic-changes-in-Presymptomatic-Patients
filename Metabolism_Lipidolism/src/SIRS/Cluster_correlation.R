# Title: SIRS correlation cluster
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

### Loading required package ---------------------------------------------------
library(vegan)
library(ggplot2)
library(reshape)
library(dplyr)
library(vegan)
library(ggpubr)
library(patchwork)
library(rstatix)
library(numbers)
library(tidyverse)
library(ggprism)
##  test  ----------------------------------------------------------------
library(vegan)
Sepsis_spearman_mat <- readRDS("Sepsis_spearman.rds")
Sepsis_r <- Sepsis_spearman_mat %>%
  dplyr::select(variable1,variable2,estimate) %>%
  pivot_wider(names_from = variable2,values_from = estimate) %>%
  column_to_rownames(var = "variable1") %>%
  as.matrix()
diag(Sepsis_r) <- 1
Sirsplus_spearman_mat <-readRDS("Sirsplus_spearman.rds")
Sirsplus_r <- Sirsplus_spearman_mat %>%
  dplyr::select(variable1,variable2,estimate) %>%
  pivot_wider(names_from = variable2,values_from = estimate) %>%
  column_to_rownames(var = "variable1") %>%
  as.matrix()
diag(Sirsplus_r) <- 1
df <- read_tsv("correlation_spearman.tsv")

Excel_data_list <- readRDS("Compound_data.rds")

Met_ID_info_hmdb <- Excel_data_list$Metabolite_data_df$compound_info %>%
  dplyr::select(CompoundID,Name) %>%
  dplyr::rename(ID = CompoundID) %>%
  mutate(ID = paste0("Met_",ID))

Lip_ID_info <- Excel_data_list$Lipid_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(LipName = Name) %>%
  dplyr::rename(ID = CompoundID) %>%
  mutate(ID = paste0("Lip_",ID))

TabS3 <- df %>%
  dplyr::rename(ID = name) %>%
  left_join(Met_ID_info_hmdb) %>%
  left_join(Lip_ID_info) %>%
  mutate(Name = if_else(is.na(Name),LipName,Name)) %>%
  dplyr::select(ID,Name,group)

library(xlsx)
write.xlsx(TabS3, "Supplementary_Tab3.xlsx",
           sheetName="Table S3",row.names = T,append=TRUE)


compound_list <- list()
for (i in c(1:5)) {
  tmp_df <- df %>%
    filter(group %in% i) %>%
    .$name 
  Sepsis_tmp <- Sepsis_r[tmp_df,tmp_df]
  Sirsplus_tmp <- Sirsplus_r[tmp_df,tmp_df]
  result_list <- list(Sepsis = Sepsis_tmp,
                      Sirsplus = Sirsplus_tmp
                      )
  compound_list[[i]] <- result_list
}
CLuster1_Sepsis.dist <- as.dist(compound_list[[1]]$Sepsis)
CLuster1_Sirsplus.dist <- as.dist(compound_list[[1]]$Sirsplus)
cor.test(CLuster1_Sepsis.dist,CLuster1_Sirsplus.dist,method = "spearman")

CLuster2_Sepsis.dist <- as.dist(compound_list[[2]]$Sepsis)
CLuster2_Sirsplus.dist <- as.dist(compound_list[[2]]$Sirsplus)
cor.test(CLuster2_Sepsis.dist,CLuster2_Sirsplus.dist,method = "spearman")

CLuster3_Sepsis.dist <- as.dist(compound_list[[3]]$Sepsis)
CLuster3_Sirsplus.dist <- as.dist(compound_list[[3]]$Sirsplus)
cor.test(CLuster3_Sepsis.dist,CLuster3_Sirsplus.dist,method = "spearman")


CLuster4_Sepsis.dist <- as.dist(compound_list[[4]]$Sepsis)
CLuster4_Sirsplus.dist <- as.dist(compound_list[[4]]$Sirsplus)
cor.test(CLuster4_Sepsis.dist,CLuster4_Sirsplus.dist,method = "spearman")


CLuster5_Sepsis.dist <- as.dist(compound_list[[5]]$Sepsis)
CLuster5_Sirsplus.dist <- as.dist(compound_list[[5]]$Sirsplus)
cor.test(CLuster5_Sepsis.dist,CLuster5_Sirsplus.dist,method = "spearman")


