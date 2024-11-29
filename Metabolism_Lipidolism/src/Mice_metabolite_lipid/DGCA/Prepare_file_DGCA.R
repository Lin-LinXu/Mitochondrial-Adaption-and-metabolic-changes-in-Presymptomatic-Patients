# Title: Data for differential correlation network - mice
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
source("functions.R")
easypackages::libraries("tidyverse", "RColorBrewer","ggpubr","rstatix",
                        "ComplexHeatmap","ppcor","patchwork","rio","missForest")
set.seed(123)
## Load data -------------------------------------------------------------------
Mice_file <- readRDS("tidy_data_new_filtered.rds")
match_lst <- readRDS("Human_signature_in_mice.rds")
Mice_meta <- import_list("RQ00754-metadata.xlsx")
Mice_meta_df <- Mice_meta$Samples %>%
  slice(-(1:4)) %>%
  janitor::row_to_names(1) %>%
  dplyr::select(`MS-Omics Sample ID`,`Treatment/ study group/ Class 1`,`Class 2`,`Class 3`) %>%
  dplyr::rename(ID = `MS-Omics Sample ID`,
                Class1 = `Treatment/ study group/ Class 1`,
                Class2 = `Class 2`,
                Class3 = `Class 3`) %>%
  mutate(Class2 = str_replace_all(Class2,c("liver \\(right lobe\\) A" = "Liver",
                                           "right kidney" = "Kidney",
                                           "spleen" = "Spleen",
                                           "WAT \\(white adipose tissue\\) A" = "WAT",
                                           "heart" = "Heart",
                                           "Serum Lipidomics" = "Serum",
                                           "Serum Metabolomics" = "Serum"
  )
  )
  ) %>%
  mutate(Class3 = substr(Class3, 1, 12)) %>%
  dplyr::select(ID,Class3)
##* Significant human signatures ----------------------------------------------- 
Ord_sig <- readRDS("Ordinal_reg_BH.rds")

sign_lip_direction <- Ord_sig$Lip_Ord$clm_df %>% 
  filter(metabolite %in% Ord_sig$Lip_Ord$re_df_025$name) %>%
  dplyr::select(metabolite,statistic,p.value,q.value) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(statistic < 0,"NonSepsis","Sepsis"))

sign_met_direction <- Ord_sig$Met_Ord$clm_df %>% 
  filter(metabolite %in% Ord_sig$Met_Ord$re_df_025$name) %>%
  dplyr::select(metabolite,statistic,p.value,q.value) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(statistic < 0,"NonSepsis","Sepsis"))

Metabolite_class_info <- read_tsv("metabolite_class_info.tsv")
Lipid_class_info <- read_tsv("Lipid_class_info.tsv")

sign_met_df <- sign_met_direction %>%
  left_join(Metabolite_class_info %>% 
              dplyr::rename(metabolite = Name) %>%
              dplyr::select(metabolite,sub_class) %>%
              dplyr::rename(Compound_class = sub_class)
  )
sign_lip_df <- sign_lip_direction %>%
  left_join(Lipid_class_info %>% 
              dplyr::rename(metabolite = Name) %>%
              dplyr::select(metabolite,Prefix) %>%
              dplyr::rename(Compound_class = Prefix)
  )

class_level <- c("Standard amino acids","Amino acids, peptides, and analogues",
                 "Ureas","Hybrid peptides","Quinoline carboxylic acids",
                 "Alpha hydroxy acids and derivatives","Imidazoles","Carbonyl compounds",
                 "lysoPC","lysoPE","TG","PC","SM","CerP","PI","PS","MG","DG","PA",
                 "plasmenyl-PE","plasmenyl-PC","DGDG")

Sign_all_signatures <- bind_rows(sign_met_df,sign_lip_df) %>%
  mutate(Compound_class = factor(Compound_class,levels = class_level))

Met_info_df <- match_lst$Human_sig_met %>%
  dplyr::select(Name,HMDB)
## Mice Statistics -------------------------------------------------------------
Tidy_mice_met <- function(df_list,match_file,Met_info){
  
  clean_met <- function(data_lst,match_df,Met_info){
    
    df_tmp <- data_lst$data_tidy_reduced_df %>%
      filter(Compound_ID %in% match_df$CompoundID) %>%
      left_join(match_df %>%
                  dplyr::rename(Compound_ID = CompoundID) %>%
                  dplyr::select(Compound_ID,HMDB)
      ) %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,HMDB,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
      )
      ) %>%
      dplyr::select(ID,Group,HMDB,value) %>%
      left_join(Met_info) %>%
      ungroup() %>%
      dplyr::select(ID,Group,Name,value)
    return(df_tmp)
  }
  
  Heart_df <- clean_met(df_list$Heart_metabolite,match_file[[2]][[1]],Met_info) %>%
    mutate(Organ = "Heart")
  Kidney_df <- clean_met(df_list$Kidney_metabolite,match_file[[2]][[2]],Met_info) %>%
    mutate(Organ = "Kidney")
  Liver_df <- clean_met(df_list$Liver_metabolite,match_file[[2]][[3]],Met_info) %>%
    mutate(Organ = "Liver")
  Spleen_df <- clean_met(df_list$Spleen_metabolite,match_file[[2]][[4]],Met_info) %>%
    mutate(Organ = "Spleen")
  WAT_df <- clean_met(df_list$Wat_metabolite,match_file[[2]][[5]],Met_info) %>%
    mutate(Organ = "WAT")
  Serum_df <- clean_met(df_list$Serum_metabolite,match_file[[2]][[6]],Met_info) %>%
    mutate(Organ = "Serum")
  
  re_df <- list(
    Heart = Heart_df,
    Kidney = Kidney_df,
    Liver = Liver_df,
    Spleen = Spleen_df,
    WAT = WAT_df,
    Serum = Serum_df
  )
  
  return(re_df)
}
Tidy_mice_lip <- function(df_list,match_file){
  
  clean_lip <- function(data_lst,match_df){
    
    df_tmp <- data_lst$data_tidy_PQN_reduced %>%
      filter(Compound_ID %in% match_df$CompoundID) %>%
      left_join(match_df %>%
                  dplyr::rename(Compound_ID = CompoundID) %>%
                  dplyr::select(Compound_ID,Short_name)
      ) %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,Short_name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
      )
      ) %>%
      group_by(ID,Group,Short_name) %>%
      summarise_at("value", mean, na.rm = TRUE) %>%
      dplyr::rename(Name = Short_name)
    return(df_tmp)
  }
  
  Heart_df <- clean_lip(df_list$Heart_lipid,match_file[[2]][[1]]) %>%
    mutate(Organ = "Heart")
  Kidney_df <- clean_lip(df_list$Kidney_lipid,match_file[[2]][[2]]) %>%
    mutate(Organ = "Kidney")
  Liver_df <- clean_lip(df_list$Liver_lipid,match_file[[2]][[3]]) %>%
    mutate(Organ = "Liver")
  Spleen_df <- clean_lip(df_list$Spleen_lipid,match_file[[2]][[4]]) %>%
    mutate(Organ = "Spleen")
  WAT_df <- clean_lip(df_list$Wat_lipid,match_file[[2]][[5]]) %>%
    mutate(Organ = "WAT")
  Serum_df <- clean_lip(df_list$Serum_lipid,match_file[[2]][[6]]) %>%
    mutate(Organ = "Serum")
  
  re_df <- list(
    Heart = Heart_df,
    Kidney = Kidney_df,
    Liver = Liver_df,
    Spleen = Spleen_df,
    WAT = WAT_df,
    Serum = Serum_df
  )
  
  return(re_df)
}

tidy_mice_met_lst <- Tidy_mice_met(Mice_file$filtered_met,match_lst$matched_met,Met_info_df)
tidy_mice_lip_lst <- Tidy_mice_lip(Mice_file$filtered_lip,match_lst$matched_lip)

## Prepare long df -------------------------------------------------------------
Combined_df <- bind_rows(
  bind_rows(tidy_mice_met_lst) %>% mutate(type = "Metabolite"),
  bind_rows(tidy_mice_lip_lst) %>% mutate(type = "Lipid")
) %>%
  left_join(Mice_meta_df) %>%
  mutate(Feature = paste0(Organ,"_",Name)) %>%
  dplyr::select(Class3,Feature,value) %>%
  pivot_wider(names_from = Feature) %>%
  column_to_rownames("Class3")
#set.seed(123)
#data_imput_list <- missForest(as.matrix(Combined_df),verbose = T)
#saveRDS(data_imput_list,file = "Mice_DGCA_impute.rds")
data_imput_list <- readRDS("Mice_DGCA_impute.rds")
data_imput_df <- data_imput_list$ximp %>%
  as_tibble() %>%
  `rownames<-`(rownames(Combined_df))

Meta_df <- bind_rows(
  bind_rows(tidy_mice_met_lst) %>% mutate(type = "Metabolite"),
  bind_rows(tidy_mice_lip_lst) %>% mutate(type = "Lipid")
) %>%
  left_join(Mice_meta_df) %>%
  dplyr::select(Class3,Group) %>%
  distinct() %>%
  dplyr::rename(Mice = Class3)

Mice_meta_DGCA <- data.frame(Mice = rownames(data_imput_df)) %>%
  left_join(Meta_df)

DGCA_prepare_lst <- list(
  data = data_imput_df,
  meta = Mice_meta_DGCA
)
saveRDS(DGCA_prepare_lst,"Mice_DGCA_prepared.rds")
