# Title: Correlation between metabolites/lipids and clinical indices
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Load library ================================================================
source("functions.R")
easypackages::libraries("tidyverse", "RColorBrewer","ggpubr","rstatix",
                        "ComplexHeatmap","ppcor","patchwork")
set.seed(123)
## Load data ===================================================================
Excel_data_list <- readRDS("Compound_data.rds")
Data_lst <- readRDS("human_data.rds")
Met_file <- Data_lst$Met_file
Lip_file <- Data_lst$Lip_file

Data_time_lst <- readRDS("human_data_time.rds")
Met_file_time <- Data_time_lst$Met_file
Lip_file_time <- Data_time_lst$Lip_file

Data_raw_lst <- readRDS("human_data_time_raw.rds")
Met_file_raw <- Data_raw_lst$Met_file
Lip_file_raw <- Data_raw_lst$Lip_file

Met_ID_info <- Excel_data_list$Metabolite_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID)
Lip_ID_info <- Excel_data_list$Lipid_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID) %>%
  mutate(Name = str_replace_all(Name, c("lysoPC" = "LPC", "lysoPE" = "LPE")))
##* Significant human signatures ----------------------------------------------- 
Ord_sig <- readRDS("Ordinal_reg_BH.rds")

sign_lip_direction <- Ord_sig$Lip_Ord$clm_df %>% 
  filter(metabolite %in% Ord_sig$Lip_Ord$re_df_025$name) %>%
  dplyr::select(metabolite,estimate) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(estimate < 0,"NonSepsis","Sepsis"))

sign_met_direction <- Ord_sig$Met_Ord$clm_df %>% 
  filter(metabolite %in% Ord_sig$Met_Ord$re_df_025$name) %>%
  dplyr::select(metabolite,estimate) %>%
  mutate(Signature = "Ord") %>%
  mutate(high_abundance = if_else(estimate < 0,"NonSepsis","Sepsis"))

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
  mutate(Compound_class = factor(Compound_class,levels = class_level)) %>%
  mutate(Compound_class = str_replace_all(Compound_class, c("lysoPC" = "LPC", "lysoPE" = "LPE")))
##* Tidy human data ------------------------------------------------------------------- 
Met_sign <- Met_file %>%
  filter(metabolite %in% sign_met_direction$metabolite) %>%
  left_join(Met_ID_info) %>%
  left_join(sign_met_direction)

Lip_sign <- Lip_file %>%
  filter(metabolite %in% sign_lip_direction$metabolite) %>%
  left_join(Lip_ID_info) %>%
  left_join(sign_lip_direction)

meta_list <- readRDS("Meta.rds")
##* Tidy human data Timepoint --------------------------------------------------
Met_sign_time <- Met_file_time %>%
  dplyr::rename(metabolite = CompoundID) %>%
  filter(metabolite %in% sign_met_direction$metabolite) %>%
  left_join(Met_ID_info) %>%
  dplyr::select(`Sample ID`,Class2,metabolite,abundance,Name,Timepoint) %>%
  left_join(sign_met_direction)

Lip_sign_time <- Lip_file_time %>%
  dplyr::rename(metabolite = CompoundID) %>%
  filter(metabolite %in% sign_lip_direction$metabolite) %>%
  left_join(Lip_ID_info) %>%
  dplyr::select(`Sample ID`,Class2,metabolite,abundance,Name,Timepoint) %>%
  left_join(sign_lip_direction)
##* Tidy meta data--------------------------------------------------------------
df_all_mean <- meta_list$df_all_mean
df_sepsis_aftTime <- meta_list$df_sepsis_aftTime

df_sepsis_bothTime <- meta_list$df_sepsis_bothTime %>%
  pivot_longer(6:length(.),names_to = "parameter",values_to = "value") %>%
  mutate(value = log2(value + 1)) %>%
  pivot_wider(names_from = Timepoint) %>%
  mutate_at(6:8, ~. - Pre) %>%
  pivot_longer(c(6:9), names_to = "Timepoint") %>%
  filter(!Timepoint %in% c("Pre","-2"))
df_sepsis_bothTime_mean <- df_sepsis_bothTime %>%
  group_by(ID,Age,Sex,parameter) %>%
  summarise_at("value", mean, na.rm = TRUE) %>%
  pivot_wider(names_from = parameter)
df_sepsis_bothTime_mean_direction <- df_sepsis_bothTime %>%
  dplyr::select(parameter) %>%
  distinct() %>%
  mutate(Sign = c(
    "neutral",
    "Sepsis",
    "NonSepsis",
    "NonSepsis",
    "NonSepsis",
    "Sepsis",
    "Sepsis",
    "Sepsis",
    "Sepsis",
    "neutral",
    "NonSepsis",
    "Sepsis",
    "Sepsis",
    "Sepsis",
    "Sepsis"
  ))

## Dataframe for correlation ---------------------------------------------------
Sepsis_Met_normalized_mean <- Met_sign_time %>%
  filter(Class2 %in% "Sepsis") %>%
  filter(Timepoint %in% c("-1","-3")) %>%
  group_by(`Sample ID`,Class2,metabolite,Name,
           Signature,high_abundance) %>%
  summarise_at("abundance", mean, na.rm = TRUE) %>%
  left_join(df_sepsis_bothTime_mean %>% dplyr::rename(`Sample ID` = ID)) %>%
  pivot_longer(c(10:24), names_to = "meta_parameter", values_to = "meta_value") %>%
  filter(Name %in% c("Serine","Glutamine","Methionine","Tryptophan",
                     "4-Acetamidobutanoic acid","N-Acetyl-alanine",
                     "Aminoadipic acid",
                     "Homocitrulline"))
Sepsis_Lip_normalized_mean <- Lip_sign_time %>%
  filter(Class2 %in% "Sepsis") %>%
  filter(Timepoint %in% c("-1","-3")) %>%
  group_by(`Sample ID`,Class2,metabolite,Name,
           Signature,high_abundance) %>%
  summarise_at("abundance", mean, na.rm = TRUE) %>%
  left_join(df_sepsis_bothTime_mean %>% dplyr::rename(`Sample ID` = ID)) %>%
  pivot_longer(c(10:24), names_to = "meta_parameter", values_to = "meta_value") %>%
  filter(Signature %in% "Ord") %>%
  left_join(Sign_all_signatures %>% dplyr::select(metabolite,Compound_class)) %>% 
  filter(Compound_class %in% c("TG","LPC","LPE","PC")) %>%
  mutate(Compound_class = factor(Compound_class,levels = c("TG","LPC","LPE","PC"))) %>%
  arrange(Compound_class) %>%
  dplyr::select(.,-Compound_class) 
# Main part --------------------------------------------------------------------
## calculate correlation 1 (normalized mean) -----------------------------------
Cor_met_cli_df <- Cal_cor(Sepsis_Met_normalized_mean,Met_ID_info)
Cor_lip_cli_df <- Cal_cor(Sepsis_Lip_normalized_mean,Lip_ID_info)

Re_cor_sepsis_normalized <- list(
  Cor_met_cli_df = Cor_met_cli_df,
  Cor_lip_cli_df = Cor_lip_cli_df
)
cor_lst <-Organize_cor_re(Re_cor_sepsis_normalized)
class_level <- c("Amino acids\nand derivatives","TG","LPC","LPE","PC")
Sign_all_signatures_df <- Sign_all_signatures %>%
  mutate(Compound_class = str_replace_all(Compound_class,c("Standard amino acids" = "Amino acids and derivatives",
                                                           "Amino acids, peptides, and analogues" = "Amino acids and derivatives")))
p <- Draw_cor_heatmap_reduced(cor_lst,df_sepsis_bothTime_mean_direction,Sign_all_signatures_df,class_level)

pdf("Cor_sepsis.pdf",
    width = 32, height = 10)
draw(p, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()


