# Title: Differential correlation network - mice
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

# Load library --------------------------------------------------------------------
library(tidyverse)
# Read data --------------------------------------------------------------------
Mice_data <- readRDS("Mice_DGCA_prepared.rds")
Human_sig_group <- readRDS("Human_signatures_group.rds")
# Tidy data --------------------------------------------------------------------
data_mat <- Mice_data$data %>%
  t()
data_mat <- Mice_data$data %>%
  select(-starts_with("Serum")) %>%
  t()
Met_design_mat <- Mice_data$meta$Group
Met_design_mat <- dummies::dummy(Met_design_mat) %>%
  `colnames<-`(c("CLP", "Sham"))
# Sign SOFA --------------------------------------------------------------------
Sign_SOFA <- readRDS("Cor_SOFA.rds")
library(DGCA)
set.seed(123)
directory_path <- "/DGCA/R/"
r_files <- list.files(directory_path, pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  source(file)
}
ddcor_Lip1000 = ddcorAll(inputMat = data_mat, design = Met_design_mat,
                         compare = c("Sham", "CLP"),
                         adjust = "perm", heatmapPlot = FALSE, 
                         nPerms = 100, 
                         corrType = "pearson")

# DGCA -------------------------------------------------------------------------
ddcor_Lip_pearson_p_dif <- ddcor_Lip1000 %>% 
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
ddcor_Lip_pearson_p_dif$Classes %>% as.factor() %>% summary()
# MEGENA cluster ---------------------------------------------------------------
set.seed(123)
ddcor_Lip1000_p_dif_megena_res_residual = ddMEGENA(ddcor_Lip_pearson_p_dif, adjusted = FALSE, evalCompactness = TRUE,pval_gene_thresh = 1,nPerm = 1000)
ddcor_Lip1000_p_dif_tmp <- ddcor_Lip_pearson_p_dif %>% filter(empPVals < 0.05)
Modules_tab <- ddcor_Lip1000_p_dif_megena_res_residual$modules %>%
  rename(Lipid = Genes)
sum_tab <- ddcor_Lip1000_p_dif_tmp$Classes %>% as.factor() %>% summary() %>% unclass()
sum_tab_re <- data.frame(x=matrix(sum_tab),row.names = names(sum_tab))


ddcor_Lip1000_module_tab <- ddcor_Lip1000_p_dif_megena_res_residual$summary$module.table$module.hub %>%
  str_split(.,"\\(\\d+\\),|\\(\\d+\\)$") %>%
  unlist()
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != ""])
ddcor_Lip1000_module_tab <- unique(ddcor_Lip1000_module_tab[ddcor_Lip1000_module_tab != "()"])
ddcor_Lip1000_module_tab_name <- str_remove_all(ddcor_Lip1000_module_tab, ".*?_") %>% unique()
write_tsv(ddcor_Lip_pearson_p_dif,file = "ddcor_pearson_mice.tsv")

## Class analyses --------------------------------------------------------------
all_group <- data.frame(Name = rownames(data_mat)) %>%
  mutate(tmp = Name) %>%
  separate(tmp, into = c("Organ","Compound"),"_") %>%
  mutate(sign = if_else(Compound %in% Sign_SOFA$Name,1,0)) %>%
  mutate(hub = if_else(Name %in% ddcor_Lip1000_module_tab,1,0)) %>%
  left_join(Human_sig_group %>% dplyr::rename(Compound = Name)) %>%
  dplyr::rename(group1 = Group)

write_tsv(all_group,file = "DGCA_mice_nodes.tsv")
