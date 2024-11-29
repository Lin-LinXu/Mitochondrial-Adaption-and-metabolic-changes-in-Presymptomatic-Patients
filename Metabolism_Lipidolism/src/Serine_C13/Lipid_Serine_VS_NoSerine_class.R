# Title: Lipid Serine+ VS Serine-
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
Tidy_mice_met_CLP <- function(df_list){

  clean_met <- function(data_lst){

    df_tmp <- data_lst$data_tidy_reduced_df %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","CLP (24h) / 12C-Serin")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "Serine-",
                                             "CLP \\(24h\\) / 12C-Serin" = "Serine+")
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
Tidy_mice_lip_CLP <- function(df_list){

  clean_lip <- function(data_lst){

    df_tmp <- data_lst$data_tidy_PQN_reduced %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","CLP (24h) / 12C-Serin")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "Serine-",
                                             "CLP \\(24h\\) / 12C-Serin" = "Serine+")
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
Tidy_mice_met_Sham <- function(df_list){

  clean_met <- function(data_lst){

    df_tmp <- data_lst$data_tidy_reduced_df %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("Sham / Vehicle","Sham / 12C-Serin")) %>%
      mutate(Group = str_replace_all(Group,c("Sham / Vehicle" = "Serine-",
                                             "Sham / 12C-Serin" = "Serine+")
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
Tidy_mice_lip_Sham <- function(df_list){

  clean_lip <- function(data_lst){

    df_tmp <- data_lst$data_tidy_PQN_reduced %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("Sham / Vehicle","Sham / 12C-Serin")) %>%
      mutate(Group = str_replace_all(Group,c("Sham / Vehicle" = "Serine-",
                                             "Sham / 12C-Serin" = "Serine+")
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

CLP_mice_met_lst <- Tidy_mice_met_CLP(Mice_file$filtered_met)
CLP_mice_lip_lst <- Tidy_mice_lip_CLP(Mice_file$filtered_lip)
Sham_mice_met_lst <- Tidy_mice_met_Sham(Mice_file$filtered_met)
Sham_mice_lip_lst <- Tidy_mice_lip_Sham(Mice_file$filtered_lip)
## T test ----------------------------------------------------------------------
T_test <- function(df_list)
{

  T_test <- function(df){

    t_df <- df %>%
      group_by(Name) %>%
      mutate(Group = factor(Group,levels = c("Serine+","Serine-"))) %>%
      t_test(value ~ Group) %>%
      mutate(Storey_p = cp4p::adjust.p(p,pi0.method="st.boot")$adjp$adjusted.p) %>%
      adjust_pvalue(method = "BH")
    
    FC_df <- df %>%
      group_by(Name,Group) %>%
      summarise_at("value",mean) %>%
      pivot_wider(names_from = Group,values_from = value) %>%
      mutate(FC = `Serine+` - `Serine-`) %>%
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

CLP_mice_lip_T_df <- T_test(CLP_mice_lip_lst)
Sham_mice_lip_T_df <- T_test(Sham_mice_lip_lst)
## Class test ----------------------------------------------------------------------
human_classlst <- c("LPC","LPE","PC","DG","PI","PS")
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

CLP_mice_lip_class_df <- Class_test(CLP_mice_lip_T_df,Mice_lip_class_lst,human_classlst)
Sham_mice_lip_class_df <- Class_test(Sham_mice_lip_T_df,Mice_lip_class_lst,human_classlst)

## Abundance table -------------------------------------------------------------
Abundance_cal <- function(df_list,class_df,human_class)
{
  
  Flt_df <- function(df,class_DF,human_classlst){
    mean_without_outliers <- function(x) {
      # Calculate Q1 (25th percentile) and Q3 (75th percentile)
      Q1 <- quantile(x, 0.25)
      Q3 <- quantile(x, 0.75)
      
      # Calculate the Interquartile Range (IQR)
      IQR <- Q3 - Q1
      
      # Define the lower and upper bounds for non-outliers
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      
      # Filter out outliers
      non_outliers <- x[x >= lower_bound & x <= upper_bound]
      
      # Calculate the mean of non-outliers
      mean(non_outliers)
    }
    
    df_new <- df %>%
      dplyr::rename(Compound_ID = Name) %>%
      left_join(class_DF) %>%
      filter(Compound_class %in% human_classlst) %>%
      group_by(Group, Compound_class) %>%
      summarise_at("value",mean) %>%
      dplyr::select(Compound_class,Group,value)
    
    return(df_new)
  }
  
  Heart_tmp <- Flt_df(df_list$Heart,class_df$Heart,human_class) %>%
    mutate(Tissue = "Heart")
  
  Kidney_tmp <- Flt_df(df_list$Kidney,class_df$Kidney,human_class) %>%
    mutate(Tissue = "Kidney")
  
  Liver_tmp <- Flt_df(df_list$Liver,class_df$Liver,human_class) %>%
    mutate(Tissue = "Liver")
  
  Spleen_tmp <- Flt_df(df_list$Spleen,class_df$Spleen,human_class) %>%
    mutate(Tissue = "Spleen")
  
  WAT_tmp <- Flt_df(df_list$WAT,class_df$WAT,human_class) %>%
    mutate(Tissue = "WAT")

  Serum_tmp <- Flt_df(df_list$Serum,class_df$Serum,human_class) %>%
    mutate(Tissue = "Serum")
  
  
  re_df <- list(
    Heart = Heart_tmp,
    Kidney = Kidney_tmp,
    Liver = Liver_tmp,
    Spleen = Spleen_tmp,
    WAT = WAT_tmp,
    Serum = Serum_tmp
  )
  
  return(re_df)
}

human_classlst <- c("LPC","LPE","PC","DG","PI","PS")

CLP_mice_lip_Abundance_df <- Abundance_cal(CLP_mice_lip_lst,Mice_lip_class_lst,human_classlst)
Sham_mice_lip_Abundance_df <- Abundance_cal(Sham_mice_lip_lst,Mice_lip_class_lst,human_classlst)

CLP_Abundance <- bind_rows(CLP_mice_lip_Abundance_df) %>% 
  mutate(Condition = "CLP")
Sham_Abundance <- bind_rows(Sham_mice_lip_Abundance_df) %>% 
  mutate(Condition = "Sham")


Total_FC <- bind_rows(CLP_Abundance,Sham_Abundance) %>%
  mutate(Type = paste0(Condition,"_",Group)) %>%
  dplyr::select(Compound_class,Tissue,Type,value) %>%
  group_by(Tissue, Compound_class) %>%
  mutate(value = scale(value)) %>%
  dplyr::select(.,-Group) %>% 
  pivot_wider(names_from = Type) %>%
  mutate(Compound_class = factor(Compound_class, levels = rev(c("LPC","LPE","PC","DG","PI","PS"))))

## Combine Plot df -------------------------------------------------------------
merge_met_lip <- function(clp_df,sham_df){
  re_df <- clp_df %>%
      dplyr::select(name,statistic,P,q) %>%
      rename_with(~ paste0("CLP_", .x), .cols = 2:ncol(.)) %>%
      left_join(sham_df %>%
                  dplyr::select(name,statistic,P,q) %>%
                  rename_with(~ paste0("Sham_", .x), .cols = 2:ncol(.)))
}
Heart_df <- merge_met_lip(CLP_mice_lip_class_df$Heart,Sham_mice_lip_class_df$Heart) %>%
  rename_with(~ paste0("Heart_", .x), .cols = 2:ncol(.))
Kidney_df <- merge_met_lip(CLP_mice_lip_class_df$Kidney,Sham_mice_lip_class_df$Kidney) %>%
  rename_with(~ paste0("Kidney_", .x), .cols = 2:ncol(.))
Liver_df <- merge_met_lip(CLP_mice_lip_class_df$Liver,Sham_mice_lip_class_df$Liver) %>%
  rename_with(~ paste0("Liver_", .x), .cols = 2:ncol(.))
Spleen_df <- merge_met_lip(CLP_mice_lip_class_df$Spleen,Sham_mice_lip_class_df$Spleen) %>%
  rename_with(~ paste0("Spleen_", .x), .cols = 2:ncol(.))
WAT_df <- merge_met_lip(CLP_mice_lip_class_df$WAT,Sham_mice_lip_class_df$WAT) %>%
  rename_with(~ paste0("WAT_", .x), .cols = 2:ncol(.))
Serum_df <- merge_met_lip(CLP_mice_lip_class_df$Serum,Sham_mice_lip_class_df$Serum) %>%
  rename_with(~ paste0("Serum_", .x), .cols = 2:ncol(.))

Plot_df <- data.frame(name = human_classlst) %>%
  left_join(Heart_df) %>%
  left_join(Kidney_df) %>%
  left_join(Liver_df) %>%
  left_join(Spleen_df) %>%
  left_join(WAT_df) %>%
  left_join(Serum_df) %>%
  mutate(name = factor(name,levels = human_classlst)) %>%
  arrange(name)

Plot_dataframe <- Plot_df %>%
  pivot_longer(-1, names_to = "parameters") %>%
  separate(parameters, into = c("Tissue", "Group", "Paramter"), sep = "_") %>%
  dplyr::rename(`Lipid Class` = name) %>%
  mutate(Paramter = str_replace_all(Paramter,c("P" = "p value","q" = "FDR")))

library(xlsx)
write.xlsx(Plot_dataframe, "Supplementary_Tab_C12.xlsx",
           sheetName="Table SX",row.names = T,append=TRUE)

FC_df <- Plot_df %>%
  dplyr::select(name, ends_with("statistic")) %>%
  rename_at(.vars = vars(ends_with("_statistic")),
            .funs = funs(sub("_statistic$", "", .))) %>%
  column_to_rownames("name") %>%
  as.matrix()

P_df <- Plot_df %>%
  dplyr::select(name,ends_with("P"))  %>%
  rename_at(.vars = vars(ends_with("_P")),
            .funs = funs(sub("_P$", "", .))) %>%
  column_to_rownames("name") %>%
  as.matrix()


FDR_df <- Plot_df %>%
  dplyr::select(name, ends_with("q")) %>%
  rename_at(.vars = vars(ends_with("_q")),
            .funs = funs(sub("_q$", "", .))) %>%
  column_to_rownames("name") %>%
  as.matrix()

FC_df[is.na(FC_df)] <- 0
P_df[is.na(P_df)] <- 1
FDR_df[is.na(FDR_df)] <- 1

P_df1 <- P_df %>%
  as.data.frame() %>%
  rownames_to_column("Compound_class") %>%
  pivot_longer(2:length(.)) %>%
  separate(name,into = c("Tissue","Group")) %>%
  mutate(value = if_else(value <= 0.05,"P","F")) %>%
  pivot_wider(names_from = Group) %>%
  dplyr::rename(P_CLP = CLP,
                P_Sham = Sham)

FDR_df1 <- FDR_df %>%
  as.data.frame() %>%
  rownames_to_column("Compound_class") %>%
  pivot_longer(2:length(.)) %>%
  separate(name,into = c("Tissue","Group")) %>%
  mutate(value = if_else(value <= 0.05,"FDR","F")) %>%
  pivot_wider(names_from = Group) %>%
  dplyr::rename(FDR_CLP = CLP,
                FDR_Sham = Sham)

Abundance_df <- Total_FC %>%
  left_join(P_df1) %>%
  left_join(FDR_df1) %>%
  mutate(Compound_class = factor(Compound_class, levels = rev(c("LPC","LPE","PC","DG","PI","PS")))) %>%
  mutate(Tissue = factor(Tissue, levels = c("Serum","WAT","Spleen","Kidney","Heart","Liver")))

p <- ggplot(Abundance_df) +
  geom_segment( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, xend=as.numeric(as.factor(Compound_class)) - 0.1, 
                    y=`CLP_Serine+`, yend=`CLP_Serine-`), color="grey") +
  geom_segment( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, xend=as.numeric(as.factor(Compound_class)) + 0.1, 
                    y=`Sham_Serine+`, yend=`Sham_Serine-`), color="grey") +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine+`,shape = 2), color="#DC0000FF", size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine-`,shape = 1), color="#DC0000FF", size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine+`), color="#3C5488FF",shape = 2, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine-`), color="#3C5488FF",shape = 1,   size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine+`, alpha= P_CLP),color="#DC0000FF",shape = 17, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine-`, alpha= P_CLP),color="#DC0000FF",shape = 16, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine+`, alpha= FDR_CLP),color="#DC0000FF",shape = 17, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) - 0.1, y=`CLP_Serine-`, alpha= FDR_CLP),color="#DC0000FF",shape = 16, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine+`, alpha= P_Sham),color="#3C5488FF",shape = 17, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine-`, alpha= P_Sham),color="#3C5488FF",shape = 16, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine+`, alpha= FDR_Sham),color="#3C5488FF",shape = 17, size=3 ) +
  geom_point( aes(x=as.numeric(as.factor(Compound_class)) + 0.1, y=`Sham_Serine-`, alpha= FDR_Sham),color="#3C5488FF",shape = 16, size=3 ) +
  scale_alpha_manual(values = c("P" = 0.3,"FDR" = 1, "F" = 0),guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = rev(c("LPC","LPE","PC","DG","PI","PS"))) +
  geom_hline(yintercept = 0, linetype="dotted", 
             color = "black", size=1) +
  facet_wrap(~Tissue, 
             scales = "free_x",
             ncol = 6) +
  coord_flip()+
  ggprism::theme_prism() +
  scale_shape_identity(name = "",
                       breaks = c(2, 1),
                       labels = c("Serine+", "Serine-"),
                       guide = "legend") +
  guides(
    color = "none",  # Remove the color legend
    shape = guide_legend(override.aes = list(color = "black"))  # Set legend shape color to black
  ) +
  xlab("") +
  ylab("Normalized mean abundance")

pdf("Serine_Compound_class.pdf",
    width = 20, height = 5)
p
dev.off()
