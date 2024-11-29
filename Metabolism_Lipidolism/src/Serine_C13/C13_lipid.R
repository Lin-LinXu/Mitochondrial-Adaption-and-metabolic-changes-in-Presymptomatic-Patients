# Title: C13 Serine lipid
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## Load Library ----------------------------------------------------------------
easypackages::libraries("tidyverse","ggprism","rstatix","ggstar","ComplexHeatmap")
jco_pal = c("CLP (24h) / Vehicle" = "#CD534CFF","CLP (24h) / 12C-Serin" = "#EFC000FF")
color_pallete <- c("CLP (24h) / Vehicle" = "#CD534CFF","Sham / Vehicle" = "#0073C2FF",
                   "CLP (24h) / 12C-Serin" = "#FF9900","Sham / 12C-Serin" = "#109618",
                   "CLP (24h) / 13C-Serin" = "#994499","Sham / 13C-Serin" = "#22AA99")
## Load Data -------------------------------------------------------------------
data <- readRDS("tidy_data_new.rds")
#* Tidy function ---------------------------------------------------------------
func_tidy_c13_lip <- function(dat_list){

  ref_ID <- dat_list$Annotation_1_df %>%
    dplyr::select(`MS-Omics ID`,`Study group`) %>%
    arrange(`MS-Omics ID`)
  
  Check_ID <- function(df,ref_ID)
  {
    check_df <- df %>%
      filter(!`Study group` %in% "QC") %>%
      dplyr::select(`MS-Omics ID`,`Study group`) %>%
      arrange(`MS-Omics ID`)
    check_re <- identical(check_df,ref_ID)
    return(check_re)
  }
  if(!Check_ID(dat_list$C13_data_df,ref_ID)){
    message("Index not same from C13 data")
  }
  
  Compound_info_df <- dat_list$C13_info_df %>%
    mutate(`Annotation level` = if_else(`Annotation level` %in% "0",
                                        "no annotation",`Annotation level`)) %>%
    rename(anno_level = `Annotation level`)
  
  Compound_data_df <- dat_list$C13_data_df %>%
    filter(!`Study group` %in% "QC") %>%
    dplyr::select(.,-`Customer ID`) %>%
    pivot_longer(3:length(.)) %>%
    filter(name %in% Compound_info_df$`Compound ID`) %>%
    pivot_wider()
  
  return_list <- list(
    Compound_info_df = Compound_info_df,
    Compound_data_df = Compound_data_df
  )
  return(return_list)
}
## Tidy data -------------------------------------------------------------------
#* Tidy Lipid -------------------------------------------------------------
Heart_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_heart_tidy_list)
Kidney_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_kidney_tidy_list)
Liver_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_liver_tidy_list)
Spleen_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_spleen_tidy_list)
Wat_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_wat_tidy_list)
Serum_lipid <- func_tidy_c13_lip(data$Lip_list$Lipid_serum_tidy_list)

Mice_C12_file <- readRDS("tidy_data_new_filtered.rds")
Tidy_lipid_class_info <- function(df)
{
  df <- df %>%
    dplyr::select(Compound_ID,Name,`Short name`) %>%
    `colnames<-`(c("Compound_ID","Name","Prefix")) %>%
    mutate(Prefix = str_remove_all(Prefix," .*")) %>%
    dplyr::select(Compound_ID,Prefix)
}
Mice_lip_class_lst <- list(
  Heart = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Heart_lipid$data_compound_info),
  Kidney = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Kidney_lipid$data_compound_info),
  Liver = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Liver_lipid$data_compound_info),
  Spleen = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Spleen_lipid$data_compound_info),
  WAT = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Wat_lipid$data_compound_info),
  Serum = Tidy_lipid_class_info(Mice_C12_file$filtered_lip$Serum_lipid$data_compound_info)
)
## Function --------------------------------------------------------------------
T_test <- function(data_list,group_vec,class_info,select_lp)
{
    
  Info_df <- data_list$Compound_info_df %>%
    rename(Compound_ID = `Compound ID`)
  
  prepare_df <- data_list$Compound_data_df %>%
    pivot_longer(3:length(.),names_to = "Compound_ID",values_to = "value") %>%
    filter(`Study group` %in% group_vec) %>%
    rename(Group = "Study group") %>%
    mutate(value = if_else(value %in% "ND","0",value)) %>%
    type_convert() %>%
    ungroup() %>%
    left_join(class_info)
    
  if (!select_lp == "all") {
    prepare_df <- prepare_df %>%
      filter(Prefix %in% select_lp)
  }else{
    prepare_df <- prepare_df %>%
      left_join(Info_df %>% dplyr::select(Compound_ID,anno_level)) %>%
      dplyr::filter(anno_level %in% c("1","2a"))
  }
  
  C13_0_vec <- prepare_df %>%
    group_by(Group,Compound_ID) %>%
    mutate(value_sum = sum(value,na.rm = T)) %>%
    filter(Group %in% group_vec[1],value_sum==0) %>%
    .$Compound_ID %>%
    unique()
  
  C12_falsePositive_vec <- prepare_df %>%
    group_by(Group,Compound_ID) %>%
    mutate(value_ave = mean(value,na.rm = T)) %>%
    filter(Group %in% group_vec[2],value_ave > 20) %>%
    .$Compound_ID %>%
    unique()

  data_t_df <- prepare_df %>%
    filter(!Compound_ID %in% C13_0_vec) %>%
    filter(!Compound_ID %in% C12_falsePositive_vec) %>%
    group_by(Compound_ID,Prefix) %>%
    mutate(value = if_else(is.na(value),0,value)) %>%
    mutate(value_variance = var(value,na.rm = T)) %>%
    filter(value_variance >= 2.6e-12) %>%
    mutate(count_val = sum(value >0,na.rm = T)) %>%
    mutate(Group = factor(Group,levels = group_vec)) %>%
    t_test(value ~ Group,alternative = "greater") %>%
    mutate(Storey_p = cp4p::adjust.p(p,pi0.method="st.boot")$adjp$adjusted.p) %>%
    adjust_pvalue(method = "BH")
  
  min_val <- prepare_df %>%
    filter(!Compound_ID %in% C13_0_vec) %>%
    filter(!Compound_ID %in% C12_falsePositive_vec) %>%
    filter(value > 0) %>%
    .$value %>%
    min()
  
  data_FC_df <- prepare_df %>%
    filter(!Compound_ID %in% C13_0_vec) %>%
    filter(!Compound_ID %in% C12_falsePositive_vec) %>%
    group_by(Compound_ID,Group) %>%
    summarise_at("value",mean) %>%
    mutate(value = log(value + min_val/2,base = 2)) %>%
    pivot_wider(names_from = Group,values_from = value) %>%
    mutate(FC = eval(as.name(group_vec[1])) - eval(as.name(group_vec[2]))) %>%
    dplyr::select(Compound_ID,FC)
  
  Class_num <- class_info %>%
    ungroup() %>%
    group_by(Prefix) %>%
    count(Prefix) %>%
    mutate(Total_num = n)
  
  data_df <- data_t_df %>%
    left_join(data_FC_df) %>%
    left_join(Info_df) %>%
    left_join(Class_num)
  
  return(data_df)
}
## Results ---------------------------------------------------------------------
group_vec <- c("CLP (24h) / 13C-Serin","CLP (24h) / 12C-Serin")
group_vec1 <- c("Sham / 13C-Serin","Sham / 12C-Serin")

T_test_Serum_clp <- T_test(Serum_lipid,group_vec,Mice_lip_class_lst$Serum,"all")
T_test_Liver_clp <- T_test(Liver_lipid,group_vec,Mice_lip_class_lst$Liver,"all")
T_test_Kidney_clp <- T_test(Kidney_lipid,group_vec,Mice_lip_class_lst$Kidney,"all")
T_test_Heart_clp <- T_test(Heart_lipid,group_vec,Mice_lip_class_lst$Heart,"all")
T_test_Spleen_clp <- T_test(Spleen_lipid,group_vec,Mice_lip_class_lst$Spleen,"all")
T_test_Wat_clp <- T_test(Wat_lipid,group_vec,Mice_lip_class_lst$WAT,"all")

T_test_Serum_sham <- T_test(Serum_lipid,group_vec1,Mice_lip_class_lst$Serum,"all")
T_test_Liver_sham <- T_test(Liver_lipid,group_vec1,Mice_lip_class_lst$Liver,"all")
T_test_Kidney_sham <- T_test(Kidney_lipid,group_vec1,Mice_lip_class_lst$Kidney,"all")
T_test_Heart_sham <- T_test(Heart_lipid,group_vec1,Mice_lip_class_lst$Heart,"all")
T_test_Spleen_sham <- T_test(Spleen_lipid,group_vec1,Mice_lip_class_lst$Spleen,"all")
T_test_Wat_sham <- T_test(Wat_lipid,group_vec1,Mice_lip_class_lst$WAT,"all")
## Enrichment ---------------------------------------------------------------------
human_classlst <- c("LPC","LPE","PC","DG","PI","PS")
Fisher_test <- function(df,class_lst){
  
  classes <- unique(class_lst)
  re_list <- lapply(classes, function(feature){
    
    subset_interest <- df %>% 
      filter(Prefix == feature & (p <= 0.05 & Storey_p <= 0.25)) %>%
      .$Compound_ID %>%
      length(.)
    
    interest_notsubset <- df %>% 
      filter(!Prefix == feature & (p <= 0.05 & Storey_p <= 0.25)) %>%
      .$Compound_ID %>%
      length(.)
      
    subset_notinterest <- df %>% 
      filter(Prefix == feature & (p > 0.05 | Storey_p > 0.25)) %>%
      .$Compound_ID %>%
      length(.)
    
    notinterest_notsubset <- df %>% 
      filter(!Prefix == feature & (p > 0.05 | Storey_p > 0.25)) %>%
      .$Compound_ID %>%
      length(.)
    
    # Create the contingency table
    contingency_table <- matrix(c(
      subset_interest,                     # DEGs in the pathway
      interest_notsubset,       # DEGs not in the pathway
      subset_notinterest,    # Non-DEGs in the pathway
      notinterest_notsubset # Non-DEGs not in the pathway
    ), nrow = 2)
    
    # Perform Fisher's Exact Test
    result <- fisher.test(contingency_table, alternative = "greater")
    
    list(name = feature, estimate = result$estimate, P = result$p.value)
  })
  
  return_df <- do.call(rbind.data.frame, re_list) %>% 
    mutate(q = p.adjust(P))
  
}

Fisher_Serum_clp <- Fisher_test(T_test_Serum_clp,human_classlst)
Fisher_Liver_clp <- Fisher_test(T_test_Liver_clp,human_classlst)
Fisher_Kidney_clp <- Fisher_test(T_test_Kidney_clp,human_classlst)
Fisher_Heart_clp <- Fisher_test(T_test_Heart_clp,human_classlst)
Fisher_Spleen_clp <- Fisher_test(T_test_Spleen_clp,human_classlst)
Fisher_Wat_clp <- Fisher_test(T_test_Wat_clp,human_classlst)

Fisher_Serum_sham <- Fisher_test(T_test_Serum_sham,human_classlst)
Fisher_Liver_sham <- Fisher_test(T_test_Liver_sham,human_classlst)
Fisher_Kidney_sham <- Fisher_test(T_test_Kidney_sham,human_classlst)
Fisher_Heart_sham <- Fisher_test(T_test_Heart_sham,human_classlst)
Fisher_Spleen_sham <- Fisher_test(T_test_Spleen_sham,human_classlst)
Fisher_Wat_sham <- Fisher_test(T_test_Wat_sham,human_classlst)

## Combine Plot df -------------------------------------------------------------
merge_met_lip <- function(clp_df,sham_df){
  re_df <- clp_df %>%
    dplyr::select(name,estimate,P,q) %>%
    rename_with(~ paste0("CLP_", .x), .cols = 2:ncol(.)) %>%
    left_join(sham_df %>%
                dplyr::select(name,estimate,P,q) %>%
                rename_with(~ paste0("Sham_", .x), .cols = 2:ncol(.)))
}
Heart_df <- merge_met_lip(Fisher_Heart_clp,Fisher_Heart_sham) %>%
  rename_with(~ paste0("Heart_", .x), .cols = 2:ncol(.))
Kidney_df <- merge_met_lip(Fisher_Kidney_clp,Fisher_Kidney_sham) %>%
  rename_with(~ paste0("Kidney_", .x), .cols = 2:ncol(.))
Liver_df <- merge_met_lip(Fisher_Liver_clp,Fisher_Liver_sham) %>%
  rename_with(~ paste0("Liver_", .x), .cols = 2:ncol(.))
Spleen_df <- merge_met_lip(Fisher_Spleen_clp,Fisher_Spleen_sham) %>%
  rename_with(~ paste0("Spleen_", .x), .cols = 2:ncol(.))
WAT_df <- merge_met_lip(Fisher_Wat_clp,Fisher_Wat_sham) %>%
  rename_with(~ paste0("WAT_", .x), .cols = 2:ncol(.))
Serum_df <- merge_met_lip(Fisher_Serum_clp,Fisher_Serum_sham) %>%
  rename_with(~ paste0("Serum_", .x), .cols = 2:ncol(.))


Tab_lst <- list(
  Fisher_Serum_sham %>% 
    mutate(Tissue = "serum", Group = "Sham"),
  Fisher_Serum_clp %>% 
    mutate(Tissue = "serum", Group = "CLP"),
  Fisher_Wat_sham %>% 
    mutate(Tissue = "white adipose tissue", Group = "Sham"),
  Fisher_Wat_clp %>% 
    mutate(Tissue = "white adipose tissue", Group = "CLP"),
  Fisher_Spleen_sham %>% 
    mutate(Tissue = "spleen", Group = "Sham"),
  Fisher_Spleen_clp %>% 
    mutate(Tissue = "spleen", Group = "CLP"),
  Fisher_Kidney_sham %>% 
    mutate(Tissue = "kidney", Group = "Sham"),
  Fisher_Kidney_clp %>% 
    mutate(Tissue = "kidney", Group = "CLP"),
  Fisher_Heart_sham %>% 
    mutate(Tissue = "heart", Group = "Sham"),
  Fisher_Heart_clp %>% 
    mutate(Tissue = "heart", Group = "CLP"),
  Fisher_Liver_sham %>% 
    mutate(Tissue = "liver", Group = "Sham"),
  Fisher_Liver_clp %>% 
    mutate(Tissue = "liver", Group = "CLP")
)

Tab_df <- bind_rows(Tab_lst) %>%
  arrange(P,q) %>%
  dplyr::rename(`Lipid Class` = name) %>%
  dplyr::rename(`Odds Ratio` = estimate) %>%
  relocate(c(Tissue,Group),.before = P)

library(xlsx)
write.xlsx(Tab_df, "Supplementary_Tab_C13.xlsx",
           sheetName="Table SX",row.names = T,append=TRUE)


