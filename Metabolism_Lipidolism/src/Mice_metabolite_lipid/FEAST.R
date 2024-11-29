# Title: FEAST analyses
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
source("functions.R")
easypackages::libraries("tidyverse", "RColorBrewer","ggpubr","rstatix",
                        "ppcor","patchwork","rio","FEAST")
set.seed(123)
## Load data -------------------------------------------------------------------
Mice_file <- readRDS("tidy_data_new_filtered.rds")
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
Mice_meta_group <- Mice_meta$Samples %>%
  slice(-(1:4)) %>%
  janitor::row_to_names(1) %>%
  dplyr::select(`MS-Omics Sample ID`,`Treatment/ study group/ Class 1`,`Class 2`,`Class 3`) %>%
  dplyr::rename(ID = `MS-Omics Sample ID`,
                Class1 = `Treatment/ study group/ Class 1`,
                Class2 = `Class 2`,
                Class3 = `Class 3`) %>%
  mutate(Class1 = str_remove_all(Class1," .*")) %>%
  mutate(Class3 = substr(Class3, 1, 12)) %>%
  dplyr::select(Class1,Class3) %>%
  distinct() %>%
  dplyr::rename(Mice = Class3,
                group = Class1)
## Tidy Data -------------------------------------------------------------------
Tidy_mice_met_raw <- function(df_list){
  
  clean_met <- function(data_lst){

    df_tmp <- data_lst$data_peak_area %>%
      pivot_longer(3:length(.), names_to = "Compound_ID", values_to = "value") %>%
      filter(!`Study group` %in% "QC") %>%
      left_join(data_lst$data_compound_info) %>%
      filter(Compound_ID %in% unique(data_lst$data_tidy_reduced_df$Compound_ID)) %>%
      filter(!is.na(InChIKey)) %>%
      dplyr::select(`MS-Omics ID`,`Study group`,Compound_ID,InChIKey,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = Compound_ID) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
      )
      ) %>%
      ungroup() %>%
      dplyr::select(ID,Group,InChIKey,value) %>%
      group_by(ID,Group,InChIKey) %>%
      summarise_at("value", mean, na.rm = TRUE) %>%
      dplyr::rename(Name = InChIKey)
    return(df_tmp)
  }
  
  Heart_df <- clean_met(df_list$Heart_metabolite) %>% mutate(Tissue = "Heart")
  Kidney_df <- clean_met(df_list$Kidney_metabolite) %>% mutate(Tissue = "Kidney")
  Liver_df <- clean_met(df_list$Liver_metabolite) %>% mutate(Tissue = "Liver")
  Spleen_df <- clean_met(df_list$Spleen_metabolite) %>% mutate(Tissue = "Spleen")
  WAT_df <- clean_met(df_list$Wat_metabolite) %>% mutate(Tissue = "WAT")
  Serum_df <- clean_met(df_list$Serum_metabolite) %>% mutate(Tissue = "Serum")
  
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
Tidy_mice_lip_raw <- function(df_list){

  clean_lip <- function(data_lst){

    df_tmp <- data_lst$data_peak_area %>%
      pivot_longer(3:length(.), names_to = "Compound_ID", values_to = "value") %>%
      filter(!`Study group` %in% "QC") %>%
      left_join(data_lst$data_compound_info) %>%
      filter(Compound_ID %in% unique(data_lst$data_tidy_PQN_reduced$Compound_ID)) %>%
      dplyr::select(`MS-Omics ID`,`Study group`,`Short name`,Name,value) %>%
      dplyr::rename(ID = `MS-Omics ID`,
                    Group = `Study group`,
                    CompoundID = `Short name`) %>%
      filter(Group %in% c("CLP (24h) / Vehicle","Sham / Vehicle")) %>%
      mutate(Group = str_replace_all(Group,c("CLP \\(24h\\) / Vehicle" = "CLP",
                                             "Sham / Vehicle" = "Sham")
      )
      ) %>%
      ungroup() %>%
      dplyr::select(ID,Group,CompoundID,value) %>%
      group_by(ID,Group,CompoundID) %>%
      summarise_at("value", mean, na.rm = TRUE) %>%
      dplyr::rename(Name = CompoundID)
    return(df_tmp)
  }
  
  Heart_df <- clean_lip(df_list$Heart_lipid) %>% mutate(Tissue = "Heart")
  Kidney_df <- clean_lip(df_list$Kidney_lipid) %>% mutate(Tissue = "Kidney")
  Liver_df <- clean_lip(df_list$Liver_lipid) %>% mutate(Tissue = "Liver")
  Spleen_df <- clean_lip(df_list$Spleen_lipid) %>% mutate(Tissue = "Spleen")
  WAT_df <- clean_lip(df_list$Wat_lipid) %>% mutate(Tissue = "WAT")
  Serum_df <- clean_lip(df_list$Serum_lipid) %>% mutate(Tissue = "Serum")
  
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

tidy_mice_met_lst <- Tidy_mice_met_raw(Mice_file$filtered_met)
tidy_mice_lip_lst <- Tidy_mice_lip_raw(Mice_file$filtered_lip)

Tidy_lst <- list(
  tidy_mice_met_lst = tidy_mice_met_lst,
  tidy_mice_lip_lst = tidy_mice_lip_lst
)
saveRDS(Tidy_lst,"Feast_tidy_data.rds")
## Functions -----------------------------------------------------------
Prepare_func <- function(data_lst,meta_df){
  
  Serum_IDs <- unique(data_lst$Serum$Name)
  data_df <- bind_rows(
    data_lst$Heart,
    data_lst$Kidney,
    data_lst$Liver,
    data_lst$Spleen,
    data_lst$WAT,
    data_lst$Serum
  ) %>%
    filter(Name %in% Serum_IDs) %>%
    left_join(meta_df) %>%
    relocate(c(Tissue,Class3),.after = Group) %>%
    pivot_wider(names_from = Name)
  
  meta_data <- data_df %>%
    filter(!Class3 %in% "Bl6J_23_0197") %>%
    dplyr::select(ID,Group,Tissue,Class3) %>%
    mutate(SourceSink = if_else(Tissue %in% "Serum","Sink","Source")) %>%
    left_join(
      data_df %>%
        filter(!Class3 %in% "Bl6J_23_0197") %>%
        ungroup() %>%
        dplyr::select(Class3) %>%
        distinct() %>%
        mutate(id = row_number()),
      by = "Class3"
    ) %>%
    mutate(Env = paste0(Tissue," ",id)) %>%
    dplyr::rename(SampleID = ID) %>%
    ungroup() %>%
    dplyr::select(SampleID,Env,SourceSink,id) %>%
    column_to_rownames("SampleID")
  
  C_mat <- data_df %>%
    filter(!Class3 %in% "Bl6J_23_0197") %>%
    ungroup() %>%
    dplyr::select(.,-c(Group,Tissue,Class3)) %>%
    column_to_rownames("ID") %>%
    as.matrix()
  
  re_list <- list(
    C_mat = C_mat,
    meta_data = meta_data
  )
  
  return(re_list)
}
Replace_NA_halfmin <- function(data_lst){

  min_value <- min(data_lst$C_mat, na.rm = TRUE)
  half_min_value <- min_value / 2
  data_lst$C_mat[is.na(data_lst$C_mat)] <- half_min_value
  
  return(data_lst)
}
FEAST_func <- function(data_lst,prefix){
  
  max_value <- max(data_lst$C_mat)
  count_matrix <- data_lst$C_mat / max_value * 1e6
  Ceiling_mat <- round(count_matrix)
  
  set.seed(123) 
  FEAST_output <- FEAST(C = Ceiling_mat, metadata = data_lst$meta_data, 
                        different_sources_flag = 1, dir_path = "/data",
                        outfile = prefix,EM_iterations = 10000)
}
Read_FEAST <- function(prefix){
  file <- paste0("/data/",prefix,"_source_contributions_matrix.txt")
  re_df <- read_tsv(file, quote = "") %>%
    rename_with(~ gsub('"', '', .)) %>%
    mutate(Name = gsub('"', '', Name)) %>%
    column_to_rownames("Name")
}
Tidy_re <- function(df,meta_df,prefix_g){
  
  re_df <- df %>%
    as.data.frame() %>%
    rownames_to_column("Name") %>%
    pivot_longer(-1) %>%
    separate(Name, into = c("Name_sink", "id1"), sep = " ") %>%
    separate(name, into = c("Name_source", "id2"), sep = " ") %>%
    filter(!is.na(value)) %>%
    mutate(Name_sink = str_replace_all(Name_sink, "_Serum", "")) %>%
    mutate(Name_source = str_replace_all(Name_source, ".*?_", "")) %>%
    dplyr::select(Name_sink,Name_source,value) %>%
    left_join(meta_df %>% dplyr::rename(Name_sink = ID)) %>%
    dplyr::rename(Mice = Class3) %>%
    mutate(Group = prefix_g)
}
Plot_stack1 <- function(df){
  
  manualcolors <- c(
    "Heart" = "#ea4335",
    "Kidney" = "#ff6d01",
    "Liver" = "#34a853",
    "Spleen" = "#4285f4",
    "WAT" = "#fbbc04",
    "Unknown" = "grey60"
  )
  
  Mice_names <- df %>%
    dplyr::select(Mice,Group) %>%
    distinct() %>%
    arrange(Group) %>%
    mutate(tmp = c(seq(1:2),seq(1:3))) %>%
    mutate(Mice_name = paste0(Group,tmp)) %>%
    dplyr::select(Mice,Mice_name)
   
  plot_df <- df %>%
    mutate(Name_source = factor(Name_source,levels = c("Heart","Kidney","Liver","Spleen",
                                                       "WAT", "Unknown")
    )
    ) %>%
    mutate(Group = factor(Group, levels = c("CLP","Sham"))) %>%
    mutate(Data = str_replace_all(Data, c("LIP" = "Lipid", "MET" = "Metabolite"))) %>%
    left_join(Mice_names)

  barchart1 <- ggplot(plot_df, aes(x = Mice_name, y = value, fill = Name_source)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ Data + Group , scales = 'free_x',  nrow = 1) +
    theme_bw() +
    scale_fill_manual(name = "Organs", values = manualcolors) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
          strip.background = element_rect(color="black", 
                                          fill="white",
                                          size= 0.3, 
                                          linetype="solid")
    )

  
  
  return(barchart1)
}
## Prepare result --------------------------------------------------------------
Met_lst <- Prepare_func(tidy_mice_met_lst,Mice_meta_df)
Lip_lst <- Prepare_func(tidy_mice_lip_lst,Mice_meta_df)

Met_halfmin <- Replace_NA_halfmin(Met_lst)
Lip_halfmin <- Replace_NA_halfmin(Lip_lst)
## Run FEAST -------------------------------------------------------------------
Met_halfmin_FEAST <- FEAST_func(Met_halfmin,"Met_halfmin")
Lip_halfmin_FEAST <- FEAST_func(Lip_halfmin,"Lip_halfmin")
## Read FEAST ------------------------------------------------------------------
Met_halfmin_FEAST_out <- Read_FEAST("Met_halfmin")
Lip_halfmin_FEAST_out <- Read_FEAST("Lip_halfmin")
## organize results ------------------------------------------------------------
Met_halfmin_FEAST_df <- Tidy_re(Met_halfmin_FEAST_out,Mice_meta_df,"MET")
Lip_halfmin_FEAST_df <- Tidy_re(Lip_halfmin_FEAST_out,Mice_meta_df,"LIP")

Halfmin_FEAST_df <- bind_rows(
  Met_halfmin_FEAST_df,
  Lip_halfmin_FEAST_df
) %>%
  left_join(Mice_meta_group) %>%
  dplyr::rename(Data = Group) %>%
  dplyr::rename(Group = group)
## Draw results ------------------------------------------------------------
Final_stack <- Plot_stack1(Halfmin_FEAST_df)
pdf("FEAST_combined.pdf",
    width = 6.5, height = 5)
Final_stack
dev.off()
