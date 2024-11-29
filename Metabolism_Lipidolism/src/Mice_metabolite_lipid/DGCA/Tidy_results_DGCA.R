# Title: Results of differential correlation network - mice
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
source("functions.R")
easypackages::libraries("tidyverse","DGCA")
set.seed(123)
## Load data -------------------------------------------------------------------
edge_file <- read_tsv("ddcor_pearson_mice.tsv")
node_file <- read_tsv("DGCA_mice_nodes.tsv")
Sign_compound <- node_file %>%
  filter(sign == 1) %>%
  .$Compound %>%
  unique()
Human_sig_group <- readRDS("Human_signatures_group.rds")
## Original file ---------------------------------------------------------------
ddcor_Lip_pearson_p_dif <- edge_file %>% 
  filter(empPVals < 0.05) %>%
  filter(!Classes %in% "0/0") %>%
  mutate(Classes = as.character(Classes)) %>%
  mutate(
    Classes = dplyr::case_when(
      Classes == "+/+" & zScoreDiff > 0 ~ "+/++",
      Classes == "+/+" & zScoreDiff < 0 ~ "++/+",
      Classes == "-/-" & zScoreDiff > 0 ~ "--/-",
      Classes == "-/-" & zScoreDiff < 0 ~ "-/--",
      TRUE                              ~ Classes
    )
  ) %>%
  filter(is.finite(zScoreDiff))

## SOFA compound ---------------------------------------------------------------
sofa_df <- node_file %>%
  filter(sign == 1)
edge_sofa_df <- ddcor_Lip_pearson_p_dif %>%
  filter(Gene1 %in% sofa_df$Name | Gene2 %in% sofa_df$Name)

sum_tab <- edge_sofa_df$Classes %>% as.factor() %>% summary() %>% unclass()
sum_tab_re <- data.frame(x=matrix(sum_tab),row.names = names(sum_tab))
sum_tab_re
set.seed(123)
ddcor_megena_res_residual = ddMEGENA(edge_sofa_df, 
                                     adjusted = FALSE, 
                                     evalCompactness = TRUE,
                                     pval_gene_thresh = 1,nPerm = 1000)

Modules_tab <- ddcor_megena_res_residual$modules %>%
  rename(Lipid = Genes)
ddcor_megena_res_residual$summary$module.table
ddcor_Lip1000_module_tab <- ddcor_megena_res_residual$summary$module.table$module.hub %>%
  str_split(.,"\\(\\d+\\),|\\(\\d+\\)$") %>%
  unlist()
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != ""])
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != "()"])
ddcor_Lip1000_module_tab_name <- str_remove_all(ddcor_Lip1000_module_tab, ".*?_") %>% unique()
ddcor_Lip1000_module_tab

sofa_node <- node_file %>%
  mutate(hub = if_else(Name %in% ddcor_Lip1000_module_tab, 1, 0)) %>%
  distinct()