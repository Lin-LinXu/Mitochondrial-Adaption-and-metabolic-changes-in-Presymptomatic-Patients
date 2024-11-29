# Title: SIRS correlation
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

# Load library -----------------------------------------------------------------
library(tidyverse)
library(qgraph)
library(janitor)
library(Hmisc)
library(ComplexHeatmap)
library(ppcor)
library(MetaboAnalystR)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(NbClust)
set.seed(123)
jco_pal = c("#8491B4FF","#DC0000FF")
value_pal = c("#5897CD","#D28636")
## Main ------------------------------------------------------------------------
# Tidy data --------------------------------------------------------------------
Excel_data_list <- readRDS("Compound_data.rds")
Data_lst <- readRDS("human_data.rds")
metabolites_hmdb_df_file <- read_tsv("Metabolite_hmdb.tsv")
metabolites_hmdb_df <- metabolites_hmdb_df_file %>%
  mutate(name = str_remove_all(name,"_r\\d$")) %>%
  distinct()
Met_file <- Data_lst$Met_file
Lip_file <- Data_lst$Lip_file
Met_ID_info_hmdb <- Excel_data_list$Metabolite_data_df$compound_info %>%
  filter(!is.na(Annotation_level)) %>%
  mutate(Molecular_Weight = as.numeric(Molecular_Weight)) %>%
  mutate(Formula = gsub(" ", "", Formula)) %>%
  left_join(metabolites_hmdb_df %>% dplyr::rename(Name = name)) %>%
  filter(CompoundID %in% unique(Met_file$metabolite))
check_hmdb <- Met_ID_info_hmdb %>%
  filter(is.na(HMDB))
Met_ID_supp <- data.frame(
  CompoundID = c("X01642","X02358","X02442","X01278","TF00252","TF00635","X00837","X00713",
           "X02442","X04004","X03297"),
  HMDB1 = c("HMDB0006548","HMDB0252111","HMDB0253028","HMDB0000897","HMDB0000508",
                  "HMDB0000452","HMDB04827","HMDB0304423","HMDB0253028","HMDB0257111",
                  "HMDB0000706")
)
Met_ID_info <- Met_ID_info_hmdb %>%
  dplyr::select(Name,CompoundID,HMDB) %>%
  left_join(Met_ID_supp) %>%
  mutate(HMDB = if_else(is.na(HMDB),HMDB1,HMDB)) %>%
  dplyr::select(.,-HMDB1) %>%
  dplyr::rename(metabolite = CompoundID) %>%
  mutate(metabolite = paste0("Met_",metabolite))
Lip_ID_info <- Excel_data_list$Lipid_data_df$compound_info %>%
  dplyr::select(Name,CompoundID) %>%
  dplyr::rename(metabolite = CompoundID) %>%
  mutate(metabolite = paste0("Lip_",metabolite))

meta_df <- read_csv2("mata_data1.csv") %>% 
  dplyr::select(.,-`Day post surgery before sepsis diagnosis or equivalent day post for comparator and SIRS Controls`) %>% 
  distinct(Classifier,Age,Gender, .keep_all = TRUE)
Clinical_df <- read_tsv("meta_data2.tsv") %>%
  mutate_at(7:length(.), as.numeric) %>%
  dplyr::select(.,-Group)
# set.seed(123)
# data_imput_list <- missForest(as.matrix(Clinical_df[,7:length(Clinical_df)]),verbose = T)
# saveRDS(data_imput_list,file = "Clinical_data_imput_list.rds")
data_imput_list <- readRDS("Clinical_data_imput_list.rds")
data_imput_df <- data_imput_list$ximp %>%
  as_tibble() %>%
  cbind(Clinical_df[,1:6],.) %>%
  pivot_longer(7:length(.),names_to = "parameter",values_to = "value") %>%
  filter(`Sample ID` %in% Met_file$`Sample ID`) 
Clinical_mat <- data_imput_df %>%
  ungroup() %>%
  dplyr::select(.,-`Surgery day difference`,-`Ethnic Origin`,-`Admission Criteria`,-Gender,-Age) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  dplyr::select(.,-Temperature)
Met_df <- Met_file %>%
  pivot_wider(names_from = metabolite, values_from = abundance)
Lip_df <- Lip_file %>%
  pivot_wider(names_from = metabolite, values_from = abundance)
same_patient <- intersect(Met_df$`Sample ID`,Lip_df$`Sample ID`)
Met_df <- Met_df %>%
  filter(`Sample ID` %in% same_patient)
Clinic_df <- Clinical_mat %>%
  filter(`Sample ID` %in% same_patient)
##Correlation Sepsis ---------------------------------------------------
## Correlation data preparation ------------------------------------------------
Met_df_Sepsis <- Met_df %>%
  filter(Class2 %in% "Sepsis") %>%
  rename_with(~ paste0("Met_", .), .cols = 3:ncol(.)) %>%
  mutate_if(is.numeric,scale)

Lip_df_Sepsis <- Lip_df %>%
  filter(Class2 %in% "Sepsis") %>%
  rename_with(~ paste0("Lip_", .), .cols = 3:ncol(.)) %>%
  mutate_if(is.numeric,scale)

Clinic_df_Sepsis <- Clinic_df %>%
  filter(`Sample ID` %in% Lip_df_Sepsis$`Sample ID`) %>%
  mutate_if(is.numeric,scale)

df_Sepsis <- cbind(Met_df_Sepsis,Lip_df_Sepsis[,-(1:2)],Clinic_df_Sepsis[,-1]) %>%
  dplyr::select(.,-Class2) %>%
  column_to_rownames("Sample ID")
df_Sepsis_tmp <- df_Sepsis %>%
  rownames_to_column()
## calculate correlation -------------------------------------------------------
corr <- rcorr(as.matrix(df_Sepsis), type = 'spearman')
r <- corr$r
p <- corr$P

Cor_df <- r %>%
  as.data.frame() %>%
  rownames_to_column("variable1") %>%
  pivot_longer(2:length(.), names_to = "variable2", values_to = "estimate")

saveRDS(Cor_df,"Sepsis_spearman.rds")
## Heatmap ---------------------------------------------------------------------
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#01befe", "white", "#ff7d00"))
hm.sepsis <- Heatmap(r, rect_gp = gpar(type = "none"), column_dend_side = "bottom",name = "Sepsis",
        show_column_names = F, show_row_names = T,
        clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
        clustering_method_columns = "complete",clustering_method_rows = "complete",
        col = col_fun,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        }
        )
km.members.sepsis <- row_order(hm.sepsis)
rows_tree_sepsis <- row_dend(hm.sepsis)
name = rownames(r)[km.members.sepsis]
cluster.index <- matrix(nrow = dim(r)[1], ncol = 1)
cluster.index[c(match(name,rownames(r))),1] <- "1"
factoextra::fviz_nbclust(r, FUN = hcut, method = "wss")
rows_clusters_sepsis <- dendextend::cutree(rows_tree_sepsis,k = 5) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  `colnames<-`(c("name","group"))
index_df_sepsis <- data.frame(name = rownames(r)[km.members.sepsis]) %>%
  left_join(rows_clusters_sepsis)
write_tsv(index_df_sepsis,"correlation_spearman.tsv")


lipid_cluster_df_sepsis <- index_df_sepsis %>%
  left_join(Lip_ID_info %>% dplyr::rename(name = metabolite)) %>%
  filter(str_starts(name, "Lip_"))
clinical_cluster_df_sepsis <- index_df_sepsis %>%
  filter(name %in% colnames(Clinic_df)) %>%
  arrange(group)

##Correlation SIRS+ ---------------------------------------------------
## Correlation data preparation ------------------------------------------------
Met_df_sirsplus <- Met_df %>%
  filter(Class2 %in% "SIRS+") %>%
  rename_with(~ paste0("Met_", .), .cols = 3:ncol(.)) %>%
  mutate_if(is.numeric,scale)

Lip_df_sirsplus <- Lip_df %>%
  filter(Class2 %in% "SIRS+") %>%
  rename_with(~ paste0("Lip_", .), .cols = 3:ncol(.)) %>%
  mutate_if(is.numeric,scale)

Clinic_df_sirsplus <- Clinic_df %>%
  filter(`Sample ID` %in% Lip_df_sirsplus$`Sample ID`) %>%
  mutate_if(is.numeric,scale)

df_Sirsplus <- cbind(Met_df_sirsplus,Lip_df_sirsplus[,-(1:2)],Clinic_df_sirsplus[,-1]) %>%
  dplyr::select(.,-Class2) %>%
  column_to_rownames("Sample ID")

df_Sirsplus_tmp <- df_Sirsplus %>%
  rownames_to_column()
## calculate correlation -------------------------------------------------------
corr <- rcorr(as.matrix(df_Sirsplus), type = 'spearman')
r <- corr$r
p <- corr$P

Cor_df <- r %>%
  as.data.frame() %>%
  rownames_to_column("variable1") %>%
  pivot_longer(2:length(.), names_to = "variable2", values_to = "estimate")

saveRDS(Cor_df,"Sirsplus_spearman.rds")

## others ----------------------------------------------------------------------
col_fun = colorRamp2(c(-1, 0, 1), c("#01befe", "white", "#ff7d00"))
draw_mat_sirsplus <- r[index_df_sepsis$name,index_df_sepsis$name]
draw_mat_sirsplus[lower.tri(draw_mat_sirsplus)] <- NA
diag(draw_mat_sirsplus) <- NA
draw_mat2 <- t(draw_mat1)
draw_mat_new <- r[index_df_sepsis$name,index_df_sepsis$name]
draw_mat_new[lower.tri(draw_mat_new)] <-  draw_mat[lower.tri(draw_mat)]
diag(draw_mat2) <- NA
diag(draw_mat1) <- NA
hm2 <- Heatmap(draw_mat_new, rect_gp = gpar(type = "none"), column_dend_side = "bottom",name = "SIRS+",
                             show_column_names = F, show_row_names = F,
                             cell_fun = function(j, i, x, y, w, h, fill) {
                               if(!is.na(draw_mat_new[i,j])) {
                                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                               }
                               if(!is.na(draw_mat1[i,j])) {
                                 grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "black"))
                               }
                               if(!is.na(draw_mat2[i,j])) {
                                 grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "black"))
                               }
                             },
                             cluster_row_slices = FALSE, cluster_column_slices = F,
                             row_split = factor(index_df_sepsis$group,levels = unique(index_df_sepsis$group)),
                             column_split = factor(index_df_sepsis$group,levels = unique(index_df_sepsis$group)),
                             row_order = index_df_sepsis$name, column_order = index_df_sepsis$name,
                             col = col_fun,
                             row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), #border = TRUE,
                             column_title = NULL,
                             row_title_gp = gpar(col = "black", fontsize = 120),
               show_heatmap_legend = F
)
pdf("correlation_spearman_cluster.pdf",
    width = 80,height = 80)
draw(hm2, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

lgd = Legend(col_fun = col_fun, title = "",legend_width = unit(4, "cm"),
             grid_height = unit(1, "cm"),
             labels_gp = gpar(col = "black", font = 84),
             direction = "horizontal")
pdf("correlation_spearman_cluster_legend.pdf",
    width = 4,height = 0.6)
draw(lgd)
dev.off()
