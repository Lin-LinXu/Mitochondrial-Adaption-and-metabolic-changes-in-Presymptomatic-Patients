# Title: Differential correlation network - human
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 27, 2024
# License: MIT

# Load library --------------------------------------------------------------------
library(DGCA)
library(tidyverse)
library(MEGENA)
# Read data --------------------------------------------------------------------
Excel_data <- readRDS("Compound_data.rds")
Lip_ID_info <- Excel_data$Lipid_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID)
Data_lst <- readRDS("human_data.rds")
Lip_file <- Data_lst$Lip_file %>%
  left_join(Lip_ID_info)
Double_names <- Lip_file %>%
  filter(`Sample ID` %in% "C033") %>%
  .$Name %>%
  table() %>%
  as.data.frame() %>%
  `colnames<-`(c("Name","Freq")) %>%
  filter(Freq > 1)
Lip_info <- Lip_file %>%
  mutate(Name = if_else(Name %in% Double_names$Name,paste0(Name,"_",metabolite),Name)) %>%
  dplyr::select(metabolite,Name) %>%
  distinct()
Lip_df <- Lip_file %>%
  mutate(Name = if_else(Name %in% Double_names$Name,paste0(Name,"_",metabolite),Name)) %>%
  dplyr::select(.,-metabolite) %>%
  pivot_wider(names_from = Name, values_from = abundance)
# read LM reg files ------------------------------------------------------------
LM_re <- readRDS("Sign_logistic_reg_result.rds")
Lip_lm <- LM_re$lipid
sign_compound <- data.frame(metabolite = Lip_lm) %>%
  left_join(Lip_info) %>%
  .$Name

Lip_mat <- Lip_df %>%
  dplyr::select(.,-Class2) %>%
  column_to_rownames(var = "Sample ID")

Met_design_mat <- Lip_df$Class2 %>% 
  gsub("Sepsis","Infection",.) %>%
  gsub("Uncomplicated Infection","Infection",.) %>%
  gsub("SIRS.*","No Infection",.)
Met_design_mat <- dummies::dummy(Met_design_mat) %>%
  `colnames<-`(c("Infection", "No Infection"))
# ------------------------------------------------------------------------------
library(DGCA)
set.seed(123)
ddcor_Lip1000 = ddcorAll(inputMat = t(Lip_mat), design = Met_design_mat,
                         compare = c("No Infection", "Infection"),
                         adjust = "perm", heatmapPlot = FALSE, nPerms = 1000, corrType = "spearman")
## Tidy results ------------------------------------------------------
ddcor_Lip_spearman_p_dif <- ddcor_Lip1000 %>% 
  filter(empPVals < 0.01) %>%
  filter(!Classes %in% "0/0") %>%
  mutate(
    Classes = case_when(
      Classes == "+/+" & zScoreDiff > 0 ~ "+/++",
      Classes == "+/+" & zScoreDiff < 0 ~ "++/+",
      Classes == "-/-" & zScoreDiff > 0 ~ "--/-",
      Classes == "-/-" & zScoreDiff < 0 ~ "-/--",
      TRUE                              ~ Classes
    )
  )
ddcor_Lip_spearman_p_dif$Classes %>% as.factor() %>% summary()
write_tsv(ddcor_Lip_spearman_p_dif,file = "Human_ddcor.tsv")

set.seed(123)
ddcor_Lip1000_p_dif_megena_res_residual = ddMEGENA(ddcor_Lip_spearman_p_dif, adjusted = FALSE, evalCompactness = TRUE,pval_gene_thresh = 1,nPerm = 1000)
saveRDS(ddcor_Lip1000_p_dif_megena_res_residual,"Human_ddcor_megena_res_residual.rds")
ddcor_Lip1000_p_dif_tmp <- ddcor_Lip_spearman_p_dif %>% filter(empPVals < 0.01)
Modules_tab <- ddcor_Lip1000_p_dif_megena_res_residual$modules %>%
  rename(Lipid = Genes)
sum_tab <- ddcor_Lip1000_p_dif_tmp$Classes %>% as.factor() %>% summary() %>% unclass()
sum_tab_re <- data.frame(x=matrix(sum_tab),row.names = names(sum_tab))

ddcor_Lip1000_module_tab <- ddcor_Lip1000_p_dif_megena_res_residual$summary$module.table$module.hub %>%
  str_split(.,"\\(\\d+\\),|\\(\\d+\\)$") %>%
  unlist()
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != ""])
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != "()"])

## Class analyses --------------------------------------------------------------
Lipid_class <- read_tsv("Lipid_class_info.tsv") %>%
  dplyr::select(Name,Prefix) %>%
  dplyr::rename(metabolite = Name,
                group1 = Prefix)
all_group <- data.frame(Name = colnames(Lip_mat),
                        shape = rep("lipid",length(colnames(Lip_mat)))
                        ) %>% 
  left_join(Lip_info) %>%
  left_join(Lipid_class) %>%
  mutate(sign = if_else(Name %in% sign_compound,"sign","nosign")) %>%
  mutate(hub = if_else(Name %in% ddcor_Lip1000_module_tab,"hub","nohub")) %>%
  mutate(group1 = if_else(is.na(group1),"FA",group1)) %>%
  mutate(group1 = if_else(group1 %in% "Fatty Acids and Conjugates","FA",group1)) %>%
  mutate(group1 = if_else(group1 %in% "Glycerophosphoglycerols","PG",group1))
write_tsv(all_group,file = "DGCA_human_nodes.tsv")
