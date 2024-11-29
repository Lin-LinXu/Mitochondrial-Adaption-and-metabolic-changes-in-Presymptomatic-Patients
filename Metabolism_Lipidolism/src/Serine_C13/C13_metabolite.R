# Title: C13 Serine metabolite
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
func_tidy_c13_met <- function(dat_list){

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
  
  if(!Check_ID(dat_list$C13_data_manual_anno_1_2a,ref_ID)){
    message("Index not same from C13 data")
  }
  if(!Check_ID(dat_list$C13_data_2b_df,ref_ID)){
    message("Index not same from C13 data")
  }
  if(!Check_ID(dat_list$C13_data_3_df,ref_ID)){
    message("Index not same from C13 data")
  }
  
  Compound_manual_info_df <- dat_list$C13_info_manual_anno_1_2a %>%
    filter(!`Precursor Adduct` %in% c("[M+H]","[M-H]")) %>%
    dplyr::select(`Compound ID`,Name) %>%
    distinct()
  
  Compound_manual_ID_df <- dat_list$C13_info_manual_anno_1_2a %>%
    filter(`Precursor Adduct` %in% c("[M+H]","[M-H]")) %>%
    dplyr::select(`Compound ID`,Name) %>%
    distinct()
  
  Compound_manual_data_df <- dat_list$C13_data_manual_anno_1_2a %>%
    pivot_longer(4:length(.)) %>%
    dplyr::rename(`Compound ID` = name) %>%
    filter(`Compound ID` %in% Compound_manual_info_df$`Compound ID`) %>%
    left_join(Compound_manual_info_df) %>%
    dplyr::select(.,-`Customer ID`,-`Compound ID`) %>%
    group_by(`MS-Omics ID`,`Study group`,Name) %>%
    type_convert() %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    left_join(Compound_manual_ID_df) %>%
    mutate(Level = "manual") %>%
    dplyr::select(`MS-Omics ID`,`Study group`,Name,`Compound ID`,Level,value)
  
  Compound_2b_data_df <- dat_list$C13_data_2b_df %>%
    pivot_longer(4:length(.)) %>%
    dplyr::rename(`Compound ID` = name) %>%
    left_join(dat_list$C13_info_2b_df %>% dplyr::select(`Compound ID`,Name)) %>%
    dplyr::select(.,-`Customer ID`) %>%
    mutate(value = str_replace_all(value, c("ND" = "0"))) %>%
    type_convert() %>%
    mutate(Level = "2b") %>%
    dplyr::select(`MS-Omics ID`,`Study group`,Name,`Compound ID`,Level,value)
  
  Compound_3_data_df <- dat_list$C13_data_3_df %>%
    pivot_longer(4:length(.)) %>%
    dplyr::rename(`Compound ID` = name) %>%
    left_join(dat_list$C13_info_3_df %>% dplyr::select(`Compound ID`,Name)) %>%
    dplyr::select(.,-`Customer ID`) %>%
    mutate(value = str_replace_all(value, c("ND" = "0"))) %>%
    type_convert() %>%
    mutate(Level = "3") %>%
    dplyr::select(`MS-Omics ID`,`Study group`,Name,`Compound ID`,Level,value) 
  
  Compound_data_df <- bind_rows(list(Compound_manual_data_df,Compound_2b_data_df,Compound_3_data_df)) %>%
    filter(!`Study group` %in% "QC") 
  
  return_list <- list(
    Compound_data_df = Compound_data_df
  )
  return(return_list)
}
## Tidy data -------------------------------------------------------------------
#* Tidy Lipid -------------------------------------------------------------
Heart_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_heart_tidy_list)
Kidney_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_kidney_tidy_list)
Liver_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_liver_tidy_list)
Spleen_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_spleen_tidy_list)
Wat_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_wat_tidy_list)
Serum_metabolite <- func_tidy_c13_met(data$Semipolar_list$Semipolar_serum_tidy_list)

## Function --------------------------------------------------------------------
T_test <- function(data_list,group_vec)
{
  
  prepare_df <- data_list$Compound_data_df %>%
    filter(`Study group` %in% group_vec) %>%
    rename(Group = "Study group") %>%
    type_convert() %>%
    ungroup() %>%
    dplyr::rename(Compound_ID = `Compound ID`)
  
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
    group_by(Compound_ID,Name,Level) %>%
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
  
  data_df <- data_t_df %>%
    left_join(data_FC_df)

  return(data_df)
}
## Results ---------------------------------------------------------------------
group_vec <- c("CLP (24h) / 13C-Serin","CLP (24h) / 12C-Serin")
group_vec1 <- c("Sham / 13C-Serin","Sham / 12C-Serin")
group_vec2 <- c("CLP (24h) / 13C-Serin","Sham / 13C-Serin")
group_vec3 <- c("Sham / 13C-Serin","CLP (24h) / 13C-Serin")

T_test_Serum_clp <- T_test(Serum_metabolite,group_vec)
T_test_Liver_clp <- T_test(Liver_metabolite,group_vec)
T_test_Kidney_clp <- T_test(Kidney_metabolite,group_vec)
T_test_Heart_clp <- T_test(Heart_metabolite,group_vec)
T_test_Spleen_clp <- T_test(Spleen_metabolite,group_vec)
T_test_Wat_clp <- T_test(Wat_metabolite,group_vec)

T_test_Serum_sham <- T_test(Serum_metabolite,group_vec1)
T_test_Liver_sham <- T_test(Liver_metabolite,group_vec1)
T_test_Kidney_sham <- T_test(Kidney_metabolite,group_vec1)
T_test_Heart_sham <- T_test(Heart_metabolite,group_vec1)
T_test_Spleen_sham <- T_test(Spleen_metabolite,group_vec1)
T_test_Wat_sham <- T_test(Wat_metabolite,group_vec1)

T_test_Serum_vec2 <- T_test(Serum_metabolite,group_vec2)
T_test_Liver_vec2 <- T_test(Liver_metabolite,group_vec2)
T_test_Kidney_vec2 <- T_test(Kidney_metabolite,group_vec2)
T_test_Heart_vec2 <- T_test(Heart_metabolite,group_vec2)
T_test_Spleen_vec2 <- T_test(Spleen_metabolite,group_vec2)
T_test_Wat_vec2 <- T_test(Wat_metabolite,group_vec2)

T_test_Serum_vec3 <- T_test(Serum_metabolite,group_vec3)
T_test_Liver_vec3 <- T_test(Liver_metabolite,group_vec3)
T_test_Kidney_vec3 <- T_test(Kidney_metabolite,group_vec3)
T_test_Heart_vec3 <- T_test(Heart_metabolite,group_vec3)
T_test_Spleen_vec3 <- T_test(Spleen_metabolite,group_vec3)
T_test_Wat_vec3 <- T_test(Wat_metabolite,group_vec3)
## Combine Plot df -------------------------------------------------------------
clp_df <- bind_rows(
  T_test_Serum_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Serum", Comparison = "CLP"),
  T_test_Liver_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Liver", Comparison = "CLP"),
  T_test_Kidney_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Kidney", Comparison = "CLP"),
  T_test_Heart_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Heart", Comparison = "CLP"),
  T_test_Spleen_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Spleen", Comparison = "CLP"),
  T_test_Wat_clp %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "WAT", Comparison = "CLP")
)
sham_df <- bind_rows(
  T_test_Serum_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Serum", Comparison = "sham"),
  T_test_Liver_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Liver", Comparison = "sham"),
  T_test_Kidney_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Kidney", Comparison = "sham"),
  T_test_Heart_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Heart", Comparison = "sham"),
  T_test_Spleen_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "Spleen", Comparison = "sham"),
  T_test_Wat_sham %>% filter(p <= 0.05, Storey_p <= 0.05) %>% mutate(Tissue = "WAT", Comparison = "sham")
)

Sign_list <- list(
  clp_df = clp_df,
  sham_df = sham_df
)

vec2_df <- bind_rows(
  T_test_Serum_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Serum", Comparison = "CLP"),
  T_test_Liver_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Liver", Comparison = "CLP"),
  T_test_Kidney_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Kidney", Comparison = "CLP"),
  T_test_Heart_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Heart", Comparison = "CLP"),
  T_test_Spleen_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Spleen", Comparison = "CLP"),
  T_test_Wat_vec2 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "WAT", Comparison = "CLP")
)
vec3_df <- bind_rows(
  T_test_Serum_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Serum", Comparison = "sham"),
  T_test_Liver_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Liver", Comparison = "sham"),
  T_test_Kidney_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Kidney", Comparison = "sham"),
  T_test_Heart_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Heart", Comparison = "sham"),
  T_test_Spleen_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "Spleen", Comparison = "sham"),
  T_test_Wat_vec3 %>% filter(p <= 0.05, Storey_p <= 0.25) %>% mutate(Tissue = "WAT", Comparison = "sham")
)

Sign_list <- list(
  More_incorporated_in_clp_df = vec2_df,
  More_incorporated_in_sham_df = vec3_df
)


clp_df_tmp <- clp_df %>%
  dplyr::select(Name,Compound_ID,Tissue,p,Storey_p) %>%
  dplyr::rename(C12_p = p, C12_Storey_p = Storey_p)

vec2_df_tmp <- vec2_df %>%
  left_join(clp_df_tmp) %>%
  filter(!is.na(C12_p))

group_vec_all <- c("Sham / 12C-Serin",
                   "Sham / 13C-Serin",
                   "CLP (24h) / 12C-Serin",
                   "CLP (24h) / 13C-Serin"
                   )

Liver_metabolite_sign_df <- Liver_metabolite$Compound_data_df %>%
  dplyr::rename(Compound_ID = `Compound ID`) %>%
  filter(Compound_ID %in% vec2_df_tmp$Compound_ID) %>%
  filter(Name %in% vec2_df_tmp$Name)  %>%
  filter(!Name %in% "Allopurinol") %>%
  filter(`Study group` %in% group_vec_all) %>%
  mutate(`Study group` = str_replace_all(`Study group`, 
                                         c("CLP \\(24h\\) / 13C-Serin" = "CLP_13C",
                                           "Sham / 13C-Serin" = "sham_13C",
                                           "CLP \\(24h\\) / 12C-Serin" = "CLP_12C",
                                           "Sham / 12C-Serin" = "sham_12C"))) %>%
  mutate(`Study group` = factor(`Study group`, levels = c("sham_12C","sham_13C",
                                                          "CLP_12C","CLP_13C"))) %>%
  mutate(Name = factor(Name)) %>%
  ungroup() %>%
  mutate(value = as.numeric(value)) %>%
  as.data.frame()

Serum_metabolite_sign_df <- Serum_metabolite$Compound_data_df %>%
  dplyr::rename(Compound_ID = `Compound ID`) %>%
  filter(Compound_ID %in% vec2_df_tmp$Compound_ID) %>%
  filter(Name %in% vec2_df_tmp$Name)  %>%
  filter(`Study group` %in% group_vec_all) %>%
  mutate(`Study group` = str_replace_all(`Study group`, 
                                         c("CLP \\(24h\\) / 13C-Serin" = "CLP_13C",
                                           "Sham / 13C-Serin" = "sham_13C",
                                           "CLP \\(24h\\) / 12C-Serin" = "CLP_12C",
                                           "Sham / 12C-Serin" = "sham_12C"))) %>%
  mutate(`Study group` = factor(`Study group`, levels = c("sham_12C","sham_13C",
                                                          "CLP_12C","CLP_13C"))) %>%
  mutate(Name = factor(Name)) %>%
  ungroup() %>%
  mutate(value = as.numeric(value)) %>%
  as.data.frame()

Combined_df <- bind_rows(
  Liver_metabolite_sign_df %>% mutate(Tissue = "Liver"),
  Serum_metabolite_sign_df %>% mutate(Tissue = "Serum")
) %>%
  mutate(Name = factor(Name))

library(ggpubr)
comparisons <- combn(unique(as.character(Combined_df$`Study group`)), 2, simplify = FALSE)
comparisons <- list(
  c("CLP_13C", "sham_13C"),
  c("CLP_13C", "CLP_12C"),
  c("sham_13C", "sham_12C")
)
library(ggprism)
jco_pal = c("sham_12C" = "#3C5488FF",
            "sham_13C" = "#EFC000FF",
            "CLP_12C" = "#DC0000FF",
            "CLP_13C" = "#868686FF")
# Create boxplots with pairwise t-tests
P <- ggplot(Combined_df, aes(x = `Study group`, y = value)) +
  geom_boxplot(aes(fill = `Study group`),outlier.shape = NA) +  # Draw boxplots without outliers
  geom_jitter( width = 0.2, alpha = 0.5, color = "black") +  # Add jitter for individual points
  facet_wrap(~ Name, scales = "free", nrow = 1) +  # Create separate plots for each compound (Name)
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif",
                     hide.ns = F, p.adjust.method = "BH") +  # Adjust p-values and hide non-significant results
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = jco_pal) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "C13 incorporation percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_prism(base_size = 16) +
  theme(
    plot.title = element_text(size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16,face= "bold",  colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black",angle = 90), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"), # bold
    strip.text.x = element_text(size = 16, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 16, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 1),
    axis.line.y = element_line(color="black", size = 1),
    legend.background=element_blank(),
    legend.text = element_text(size=16, face= "bold"),
    legend.title = element_text(size=16, face= "bold")
  )

pdf("C13_metabolit.pdf",
    width = 25,height = 5)
P
dev.off()

