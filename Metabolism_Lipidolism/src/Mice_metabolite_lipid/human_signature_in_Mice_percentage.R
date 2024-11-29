# Title: Human signatures in mice
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
library(tidyverse)
library(ggvenn)
library(patchwork)
## Load data -------------------------------------------------------------------
Ord_sig <- readRDS("Ordinal_reg_BH.rds")

Excel_data_list <- readRDS("Compound_data.rds")
metabolites_hmdb_df <- read_tsv("Metabolite_hmdb.tsv")

Met_ID_info <- Excel_data_list$Metabolite_data_df$compound_info %>%
  filter(!is.na(Annotation_level)) %>%
  mutate(Molecular_Weight = as.numeric(Molecular_Weight)) %>%
  mutate(Formula = gsub(" ", "", Formula)) %>%
  left_join(metabolites_hmdb_df %>% dplyr::rename(Name = name))
Lip_ID_info <- Excel_data_list$Lipid_data_df$compound_info %>%
  filter(!is.na(Annotation_level)) %>%
  mutate(tmp = Name) %>%
  separate(tmp,into=c("Short_name","Others"),";") %>%
  mutate(Short_name = replace_na(Short_name,"NA_chr")) %>%
  mutate(Short_name = str_replace_all(Short_name,c("lysoPC" = "LPC","lysoPE" = "LPE","Ceramide" = "Cer"))) %>%
  dplyr::select(.,-Others) %>%
  mutate(Molecular_Weight = as.numeric(Molecular_Weight)) %>%
  mutate(Formula = gsub(" ", "", Formula))

Mice_file <- readRDS("tidy_data_new_filtered.rds")
## Organize human results ------------------------------------------------------------
Human_sig_met <- data.frame(
  CompoundID = unique(Ord_sig$Met_Ord$re_df_025$name) 
) %>%
  left_join(Met_ID_info) %>%
  dplyr::select(CompoundID,Name,RT,Molecular_Weight,Formula,HMDB)
Human_sig_lip <- data.frame(
  CompoundID = unique(c(Ord_sig$Lip_Ord$re_df_025$name))
) %>%
  left_join(Lip_ID_info) %>%
  dplyr::select(CompoundID,Name,RT,Molecular_Weight,Formula,Short_name)

Human_sig_profile <- bind_rows(
  data.frame(Name = Human_sig_met$Name) %>% mutate(Group = "Metabolite"),
  data.frame(Name = Human_sig_lip$Short_name) %>% 
    mutate(Group = str_remove_all(Human_sig_lip$Short_name," .*"))
)

total_lip <- Human_sig_profile %>%
  filter(!Group %in% "Metabolite") %>%
  dplyr::select(Name) %>%
  distinct()
## Clean mice info -------------------------------------------------------------------------
Clean_mice_met <- function(df_list)
{
  clean_met <- function(data_lst){
    df_info <- data_lst$data_compound_info %>%
      filter(Compound_ID %in% data_lst$data_tidy_reduced_df$Compound_ID) %>%
      dplyr::select(Compound_ID,Name,Retention_time,Molecular_weight,Formula,HMDB_ID) %>%
      dplyr::rename(CompoundID = Compound_ID,
                    RT = Retention_time,
                    Molecular_Weight = Molecular_weight,
                    HMDB = HMDB_ID)
  }
  Heart_info <- clean_met(df_list$Heart_metabolite)
  Kidney_info <- clean_met(df_list$Kidney_metabolite)
  Liver_info <- clean_met(df_list$Liver_metabolite)
  Spleen_info <- clean_met(df_list$Spleen_metabolite)
  WAT_info <- clean_met(df_list$Wat_metabolite)
  Serum_info <- clean_met(df_list$Serum_metabolite)
  re_info <- list(
    Heart = Heart_info,
    Kidney = Kidney_info,
    Liver = Liver_info,
    Spleen = Spleen_info,
    WAT = WAT_info,
    Serum = Serum_info
  )
}
Clean_mice_lip <- function(df_list)
{
  clean_met <- function(data_lst){
    df_info <- data_lst$data_compound_info %>%
      filter(Compound_ID %in% data_lst$data_tidy_PQN_reduced$Compound_ID) %>%
      dplyr::select(Compound_ID,Name,Retention_time,Molecular_weight,Formula,`Short name`) %>%
      dplyr::rename(CompoundID = Compound_ID,
                    RT = Retention_time,
                    Molecular_Weight = Molecular_weight,
                    Short_name = `Short name`)
  }
  Heart_info <- clean_met(df_list$Heart_lipid)
  Kidney_info <- clean_met(df_list$Kidney_lipid)
  Liver_info <- clean_met(df_list$Liver_lipid)
  Spleen_info <- clean_met(df_list$Spleen_lipid)
  WAT_info <- clean_met(df_list$Wat_lipid)
  Serum_info <- clean_met(df_list$Serum_lipid)
  re_info <- list(
    Heart = Heart_info,
    Kidney = Kidney_info,
    Liver = Liver_info,
    Spleen = Spleen_info,
    WAT = WAT_info,
    Serum = Serum_info
  )
}
Mice_met_info <- Clean_mice_met(Mice_file$filtered_met)
Mice_lip_info <- Clean_mice_lip(Mice_file$filtered_lip)

## Matching human and mice -----------------------------------------------------
Match_mice_met <- function(mice_lst,human_df)
{

  match_met <- function(data_df){
    df_info <- data_df %>%
      filter(HMDB %in% human_df$HMDB |
               Name %in% human_df$Name)
  }
  Heart_info <- match_met(mice_lst$Heart)
  Kidney_info <- match_met(mice_lst$Kidney)
  Liver_info <- match_met(mice_lst$Liver)
  Spleen_info <- match_met(mice_lst$Spleen)
  WAT_info <- match_met(mice_lst$WAT)
  Serum_info <- match_met(mice_lst$Serum)
  
  re_tb <- profs <- tribble(
    ~ Tissue, ~ Profile, ~ Num, 
    "Heart", Heart_info, dim(Heart_info)[1],
    "Kidney", Kidney_info, dim(Kidney_info)[1],
    "Liver", Liver_info, dim(Liver_info)[1],
    "Spleen", Spleen_info, dim(Spleen_info)[1],
    "WAT", WAT_info, dim(WAT_info)[1],
    "Serum", Serum_info, dim(Serum_info)[1]
  ) %>%
    mutate(Percentage = round(Num/dim(human_df)[1],2) * 100) %>%
    mutate(Percentage = paste0(Percentage,"%"))

}
Match_mice_lip <- function(mice_lst,human_df)
{

  match_lip <- function(data_df){
    df_info <- data_df %>%
      filter(Short_name %in% human_df$Short_name |
               Name %in% human_df$Name)
  }
  Heart_info <- match_lip(mice_lst$Heart)
  Kidney_info <- match_lip(mice_lst$Kidney)
  Liver_info <- match_lip(mice_lst$Liver)
  Spleen_info <- match_lip(mice_lst$Spleen)
  WAT_info <- match_lip(mice_lst$WAT)
  Serum_info <- match_lip(mice_lst$Serum)
  
  re_tb <- profs <- tribble(
    ~ Tissue, ~ Profile, ~ Num, 
    "Heart", Heart_info, length(unique(Heart_info$Short_name)),
    "Kidney", Kidney_info, length(unique(Kidney_info$Short_name)),
    "Liver", Liver_info, length(unique(Liver_info$Short_name)),
    "Spleen", Spleen_info, length(unique(Spleen_info$Short_name)),
    "WAT", WAT_info, length(unique(WAT_info$Short_name)),
    "Serum", Serum_info, length(unique(Serum_info$Short_name))
  ) %>%
    mutate(Percentage = round(Num/length(unique(human_df$Short_name)),2) * 100) %>%
    mutate(Percentage = paste0(Percentage,"%"))
  
}

## Matching results ------------------------------------------------------------
matched_met <- Match_mice_met(Mice_met_info,Human_sig_met)
matched_lip <- Match_mice_lip(Mice_lip_info,Human_sig_lip)

met_all <- bind_rows(matched_met$Profile) %>%
  dplyr::select(HMDB) %>%
  distinct()

met_all_df <- data.frame(
  Tissue = "All",
  Num = dim(met_all)[1],
  Percentage = round(dim(met_all)[1]/15,2) * 100 
) %>%
  mutate(Percentage = paste0(Percentage,"%"))

lip_all <- bind_rows(matched_lip$Profile) %>%
  dplyr::select(Short_name) %>%
  distinct()

lip_all_df <- data.frame(
  Tissue = "All",
  Num = dim(lip_all)[1],
  Percentage = round(dim(lip_all)[1]/93,2) * 100 
) %>%
  mutate(Percentage = paste0(Percentage,"%"))

re_match_lst <- list(
  matched_met = matched_met,
  matched_lip = matched_lip,
  Human_sig_met = Human_sig_met,
  Human_sig_lip = Human_sig_lip
)
saveRDS(re_match_lst,"Human_signature_in_mice.rds")

matched_met <- Match_mice_met(Mice_met_info,Human_sig_met) %>%
  dplyr::select(.,-Profile) %>%
  arrange(Num) %>%
  mutate(Tissue = factor(Tissue,levels = .$Tissue)) %>%
  rbind(met_all_df)
matched_lip <- Match_mice_lip(Mice_lip_info,Human_sig_lip) %>%
  dplyr::select(.,-Profile) %>%
  arrange(Num) %>%
  mutate(Tissue = factor(Tissue,levels = .$Tissue)) %>%
  rbind(lip_all_df)

## Final figure ----------------------------------------------------------------
Combined_matched <- bind_rows(
  matched_met %>% mutate(type = "metabolite"),
  matched_lip %>% mutate(type = "lipid")
) %>%
  mutate(Percentage_num = str_remove_all(Percentage,"%"))
library(ggprism)
P <- ggplot(Combined_matched, aes(x = Tissue, y = Num, fill = type, group = type)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("metabolite" = "#4DBBD5FF", "lipid" = "#00A087FF")) + 
  labs(title = "",
       x = "",
       y = "Number of signature compounds",
       fill = "") +
  coord_flip() +
  geom_text(aes(label = Percentage), size = 6,
            fontface = "bold",
            position = position_dodge(width = 0.8)
            ) +
  theme_classic() +
  theme_prism(base_size = 16) +
  theme(legend.position = "bottom",legend.text=element_text(size=16, face = "bold"))


P

pdf("Human_signature_in_mice.pdf",
    width = 6,height = 8)
P
dev.off()


