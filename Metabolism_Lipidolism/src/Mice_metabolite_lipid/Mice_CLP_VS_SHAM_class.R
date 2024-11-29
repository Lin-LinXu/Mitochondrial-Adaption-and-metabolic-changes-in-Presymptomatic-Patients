# Title: Mice CLP VS Sham
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
source("functions.R")
easypackages::libraries("tidyverse", "RColorBrewer","ggpubr","rstatix",
                        "ComplexHeatmap","ppcor","patchwork")
set.seed(123)
## Load data -------------------------------------------------------------------
Mice_file <- readRDS("tidy_data_new_filtered.rds")
##* Significant human signatures ----------------------------------------------- 
Ord_sig <- readRDS("Ordinal_reg_BH.rds")

sign_lip_direction <- Ord_sig$Lip_Ord$clm_df %>% 
  dplyr::select(metabolite,statistic) %>%
  dplyr::rename(estimate = statistic) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(estimate < 0,"NonSepsis","Sepsis"))

sign_met_direction <- Ord_sig$Met_Ord$clm_df %>% 
  dplyr::select(metabolite,statistic) %>%
  dplyr::rename(estimate = statistic) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(estimate < 0,"NonSepsis","Sepsis"))

Metabolite_class_info <- read_tsv("metabolite_class_info.tsv")
Lipid_class_info <- read_tsv("Lipid_class_info.tsv")

sign_met_df <- sign_met_direction %>%
  left_join(Metabolite_class_info %>% 
              dplyr::rename(metabolite = Name) %>%
              dplyr::select(metabolite,sub_class) %>%
              dplyr::rename(Compound_class = sub_class)
  ) %>%
  filter(!is.na(Compound_class))
sign_lip_df <- sign_lip_direction %>%
  left_join(Lipid_class_info %>% 
              dplyr::rename(metabolite = Name) %>%
              dplyr::select(metabolite,Prefix) %>%
              dplyr::rename(Compound_class = Prefix)
  ) %>%
  filter(!is.na(Compound_class)) %>%
  mutate(Compound_class = str_replace_all(Compound_class,
                                          c("lysoPE" = "LPE",
                                            "lysoPC" = "LPC",
                                            "Ceramide" = "Cer",
                                            "Fatty Acids and Conjugates" = "FA")))

##* Mice classes ---------------------------------------------------------------
Class_info_list <- readRDS("class_info_list.rds")

Amino_acid_list <- c("Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
                     "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
                     "Leucine","Lysine","Methionine","Phenylalanine","Proline",
                     "Serine","Threonine","Tryptophan","Tyrosine","Valine"
)

Reassign_class_met <- function(df,class_df,Amino_acid_list){
  re_df <- df %>%
    left_join(class_df) %>%
    dplyr::select(Compound_ID,Name,Subclass) %>%
    mutate(Subclass = if_else(Name %in% Amino_acid_list,"Standard amino acids",Subclass)) %>%
    filter(!is.na(Subclass)) %>%
    dplyr::rename(Compound_class = Subclass)
}

Tidy_lipid_class_info <- function(df)
{
  df <- df %>%
    dplyr::select(Compound_ID,Name,`Short name`) %>%
    `colnames<-`(c("Compound_ID","Name","Prefix")) %>%
    mutate(Prefix = str_remove_all(Prefix," .*")) %>%
    dplyr::rename(Compound_class = Prefix)
}

Mice_met_class_lst <- list(
  Heart = Reassign_class_met(Mice_file$filtered_met$Heart_metabolite$data_compound_info,
                             Class_info_list$Metabolite_class_info_list$Heart_metabolite_class_final,
                             Amino_acid_list),
  Kidney = Reassign_class_met(Mice_file$filtered_met$Kidney_metabolite$data_compound_info,
                             Class_info_list$Metabolite_class_info_list$Kidney_metabolite_class_final,
                             Amino_acid_list),
  Liver = Reassign_class_met(Mice_file$filtered_met$Liver_metabolite$data_compound_info,
                              Class_info_list$Metabolite_class_info_list$Liver_metabolite_class_final,
                              Amino_acid_list),
  Spleen = Reassign_class_met(Mice_file$filtered_met$Spleen_metabolite$data_compound_info,
                             Class_info_list$Metabolite_class_info_list$Spleen_metabolite_class_final,
                             Amino_acid_list),
  WAT = Reassign_class_met(Mice_file$filtered_met$Wat_metabolite$data_compound_info,
                              Class_info_list$Metabolite_class_info_list$Wat_metabolite_class_final,
                              Amino_acid_list),
  Serum = Reassign_class_met(Mice_file$filtered_met$Serum_metabolite$data_compound_info,
                           Class_info_list$Metabolite_class_info_list$Serum_metabolite_class_final,
                           Amino_acid_list)
)

Mice_lip_class_lst <- list(
  Heart = Tidy_lipid_class_info(Mice_file$filtered_lip$Heart_lipid$data_compound_info),
  Kidney = Tidy_lipid_class_info(Mice_file$filtered_lip$Kidney_lipid$data_compound_info),
  Liver = Tidy_lipid_class_info(Mice_file$filtered_lip$Liver_lipid$data_compound_info),
  Spleen = Tidy_lipid_class_info(Mice_file$filtered_lip$Spleen_lipid$data_compound_info),
  WAT = Tidy_lipid_class_info(Mice_file$filtered_lip$Wat_lipid$data_compound_info),
  Serum = Tidy_lipid_class_info(Mice_file$filtered_lip$Serum_lipid$data_compound_info)
)
## Mice Statistics -------------------------------------------------------------
Tidy_mice_met <- function(df_list){

  clean_met <- function(data_lst){

    df_tmp <- data_lst$data_tidy_reduced_df %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
                                     )
      ) %>%
      dplyr::select(ID,Group,CompoundID,value) %>%
      dplyr::rename(Name = CompoundID)
    return(df_tmp)
  }
  
  Heart_df <- clean_met(df_list$Heart_metabolite)
  Kidney_df <- clean_met(df_list$Kidney_metabolite)
  Liver_df <- clean_met(df_list$Liver_metabolite)
  Spleen_df <- clean_met(df_list$Spleen_metabolite)
  WAT_df <- clean_met(df_list$Wat_metabolite)
  Serum_df <- clean_met(df_list$Serum_metabolite)
  
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
Tidy_mice_lip <- function(df_list){

  clean_lip <- function(data_lst){

    df_tmp <- data_lst$data_tidy_PQN_reduced %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
      )
      ) %>%
      dplyr::select(ID,Group,CompoundID,value) %>%
      dplyr::rename(Name = CompoundID)
    return(df_tmp)
  }
  
  Heart_df <- clean_lip(df_list$Heart_lipid)
  Kidney_df <- clean_lip(df_list$Kidney_lipid)
  Liver_df <- clean_lip(df_list$Liver_lipid)
  Spleen_df <- clean_lip(df_list$Spleen_lipid)
  WAT_df <- clean_lip(df_list$Wat_lipid)
  Serum_df <- clean_lip(df_list$Serum_lipid)
  
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

tidy_mice_met_lst <- Tidy_mice_met(Mice_file$filtered_met)
tidy_mice_lip_lst <- Tidy_mice_lip(Mice_file$filtered_lip)
## T test ----------------------------------------------------------------------
T_test <- function(df_list)
{

  T_test <- function(df){

    t_df <- df %>%
      group_by(Name) %>%
      t_test(value ~ Group) %>%
      mutate(Storey_p = cp4p::adjust.p(p,pi0.method="st.boot")$adjp$adjusted.p) %>%
      adjust_pvalue(method = "BH")
    
    FC_df <- df %>%
      group_by(Name,Group) %>%
      summarise_at("value",mean) %>%
      pivot_wider(names_from = Group,values_from = value) %>%
      mutate(FC = CLP - Sham) %>%
      dplyr::select(Name,FC)
    
    com_df <- t_df %>%
      left_join(FC_df)
    
    return(com_df)
  }
  
  Heart_df <- T_test(df_list$Heart)
  Kidney_df <- T_test(df_list$Kidney)
  Liver_df <- T_test(df_list$Liver)
  Spleen_df <- T_test(df_list$Spleen)
  WAT_df <- T_test(df_list$WAT)
  Serum_df <- T_test(df_list$Serum)
  
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

mice_met_T_df <- T_test(tidy_mice_met_lst)
mice_lip_T_df <- T_test(tidy_mice_lip_lst)
## Combine preparation -------------------------------------------------------------
Class_freq_met <- table(sign_met_df$Compound_class) %>%
  as.data.frame() %>%
  filter(Freq > 5)
Class_freq_lip <- table(sign_lip_df$Compound_class) %>%
  as.data.frame() %>%
  filter(Freq > 5)

class_level <- c("Standard amino acids","Amino acids, peptides, and analogues",
                 "Carbohydrates and carbohydrate conjugates","Carbonyl compounds",
                 "Fatty acid esters", "Fatty acids and conjugates",
                 "Purines and purine derivatives",
                 "LPC","LPE","PC","TG","DG","SM","Cer","CerP","PS","PI","PA","FA",
                 "plasmenyl-PE","plasmenyl-PC")

Sign_all_signatures_class <- bind_rows(Class_freq_met,Class_freq_lip) %>%
  dplyr::rename(Compound_class = Var1) %>%
  dplyr::select(.,-Freq) %>%
  filter(Compound_class %in% class_level)

Human_estimate_list <- list(
  met_df = sign_met_df %>% filter(Compound_class %in% class_level),
  lip_df = sign_lip_df %>% filter(Compound_class %in% class_level)
)
## Class test ----------------------------------------------------------------------
Class_test_human <- function(df)
{

  T_one_VS_others <- function(df){
    classes <- unique(df$Compound_class)
    re_list <- lapply(classes, function(feature){
      x <- df %>% filter(Compound_class == feature) %>% .$estimate %>% as.numeric()
      list(name = feature, statistic = t.test(x)$statistic, P = t.test(x)$p.value)
    })
    return_df <- do.call(rbind.data.frame, re_list)
    return(return_df)
  }
  
  re_df <- T_one_VS_others(df) %>% 
    mutate(q = p.adjust(P))
  
  return(re_df)
}
Class_test <- function(df_list,class_df,human_class)
{
  
  Flt_df <- function(df,class_DF,human_classlst){
  
    df_new <- df %>%
      dplyr::rename(Compound_ID = Name) %>%
      left_join(class_DF) %>%
      filter(Compound_class %in% human_classlst) %>%
      dplyr::select(Compound_ID,statistic,Compound_class)
    
    df_freq <- table(df_new$Compound_class) %>%
      as.data.frame() %>%
      filter(Freq >= 2)

    df_new <- df_new %>%
      filter(Compound_class %in% df_freq$Var1)
    
    return(df_new)
  }
  
  T_one_VS_others <- function(df){
    classes <- unique(df$Compound_class)
    re_list <- lapply(classes, function(feature){
      x <- df %>% filter(Compound_class == feature) %>% .$statistic %>% as.numeric()
      list(name = feature, statistic = t.test(x)$statistic, P = t.test(x)$p.value)
    })
    return_df <- do.call(rbind.data.frame, re_list)
    return(return_df)
  }
  
  Heart_tmp <- Flt_df(df_list$Heart,class_df$Heart,human_class)
  Heart_df <- T_one_VS_others(Heart_tmp) %>% 
    mutate(q = p.adjust(P))
  
  Kidney_tmp <- Flt_df(df_list$Kidney,class_df$Kidney,human_class)
  Kidney_df <- T_one_VS_others(Kidney_tmp) %>% 
    mutate(q = p.adjust(P))
  
  Liver_tmp <- Flt_df(df_list$Liver,class_df$Liver,human_class)
  Liver_df <- T_one_VS_others(Liver_tmp) %>% 
    mutate(q = p.adjust(P))
  
  Spleen_tmp <- Flt_df(df_list$Spleen,class_df$Spleen,human_class)
  Spleen_df <- T_one_VS_others(Spleen_tmp) %>% 
    mutate(q = p.adjust(P))
  
  WAT_tmp <- Flt_df(df_list$WAT,class_df$WAT,human_class)
  WAT_df <- T_one_VS_others(WAT_tmp) %>% 
    mutate(q = p.adjust(P))
  
  Serum_tmp <- Flt_df(df_list$Serum,class_df$Serum,human_class)
  Serum_df <- T_one_VS_others(Serum_tmp) %>% 
    mutate(q = p.adjust(P))
  
  
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
human_class_met_df <- Class_test_human(Human_estimate_list$met_df)
human_class_lip_df <- Class_test_human(Human_estimate_list$lip_df)
mice_met_class_df <- Class_test(mice_met_T_df,Mice_met_class_lst,class_level)
mice_lip_class_df <- Class_test(mice_lip_T_df,Mice_lip_class_lst,class_level)
## Combine Plot df -------------------------------------------------------------
merge_met_lip <- function(met_df,lip_df){
  re_df <- bind_rows(
    met_df %>%
      dplyr::select(name,statistic,P,q),
    lip_df %>%
      dplyr::select(name,statistic,P,q)
  )
}
Heart_df <- merge_met_lip(mice_met_class_df$Heart,mice_lip_class_df$Heart) %>%
  rename_with(~ paste0("Heart_", .x), .cols = 2:ncol(.))
Kidney_df <- merge_met_lip(mice_met_class_df$Kidney,mice_lip_class_df$Kidney) %>%
  rename_with(~ paste0("Kidney_", .x), .cols = 2:ncol(.))
Liver_df <- merge_met_lip(mice_met_class_df$Liver,mice_lip_class_df$Liver) %>%
  rename_with(~ paste0("Liver_", .x), .cols = 2:ncol(.))
Spleen_df <- merge_met_lip(mice_met_class_df$Spleen,mice_lip_class_df$Spleen) %>%
  rename_with(~ paste0("Spleen_", .x), .cols = 2:ncol(.))
WAT_df <- merge_met_lip(mice_met_class_df$WAT,mice_lip_class_df$WAT) %>%
  rename_with(~ paste0("WAT_", .x), .cols = 2:ncol(.))
Serum_df <- merge_met_lip(mice_met_class_df$Serum,mice_lip_class_df$Serum) %>%
  rename_with(~ paste0("Serum_", .x), .cols = 2:ncol(.))

Plot_df <- bind_rows(human_class_met_df,human_class_lip_df) %>%
  rename_with(~ paste0("Human_", .x), .cols = 2:ncol(.)) %>%
  left_join(Heart_df) %>%
  left_join(Kidney_df) %>%
  left_join(Liver_df) %>%
  left_join(Spleen_df) %>%
  left_join(WAT_df) %>%
  left_join(Serum_df) %>%
  mutate(name = factor(name,levels = class_level)) %>%
  arrange(name)

saveRDS(Plot_df,"all_compounds_class_T.rds")

