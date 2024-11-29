# Title: Mice CLP VS Sham plot
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

## load library ----------------------------------------------------------------
source("functions.R")
easypackages::libraries("tidyverse", "RColorBrewer","ggpubr","rstatix",
                        "ComplexHeatmap","ppcor","patchwork")
set.seed(123)
## Prepare Heatmap ----------------------------------------------------------------
Human_Signatures_T <- readRDS("Human_Signatures_T.rds")
Class_T <- readRDS("Mice_Signatures_T.rds")

metabolite_class_level <- c("Standard amino acids","Amino acids, peptides, and analogues",
                 "Ureas","Hybrid peptides","Quinoline carboxylic acids",
                 "Alpha hydroxy acids and derivatives","Imidazoles","Carbonyl compounds")
lipid_class_level <- c("LPC","LPE","PC","PI","DG")

Metabolite_df <- Human_Signatures_T %>%
  filter(Compound_class %in% metabolite_class_level) %>%
  dplyr::select(.,-c(match_name,Signature,high_abundance)) %>%
  dplyr::select(-matches("Storey_p$")) %>%
  rename_all(~ sub("_p$", "_P", .)) %>%
  rename_all(~ sub("_p.adj$", "_q", .)) %>%
  dplyr::rename(Human_statistic = statistic,
                Human_P = `p.value`,
                Human_q = `q.value`
                ) %>%
  mutate(Compound_class = factor(Compound_class,metabolite_class_level)) %>%
  arrange(Compound_class) %>%
  filter(Compound_class %in% c("Standard amino acids","Amino acids, peptides, and analogues")) %>%
  ungroup() %>%
  dplyr::select(.,-c(match_name,Signature,Compound_class)) %>%
  mutate(Group = "amino acids and amino acid derivatives")
lipid_df <- Class_T %>%
  filter(name %in% lipid_class_level) %>%
  mutate(Group = "Lipid")

Metabolite_df <- bind_rows(slice(Metabolite_df, 3), slice(Metabolite_df, -3))
Plot_df <- bind_rows(lipid_df,Metabolite_df)

## Heatmap plot ----------------------------------------------------------------
FC_df <- Plot_df %>%
  dplyr::select(name, ends_with("statistic")) %>%
  rename_at(.vars = vars(ends_with("_statistic")),
            .funs = funs(sub("_statistic$", "", .))) %>%
  column_to_rownames("name") %>%
  relocate(Serum,.after = "Human") %>%
  rename_with(~ ifelse(. == "Human", "Human\nserum", .)) %>%
  as.matrix()

P_df <- Plot_df %>%
  dplyr::select(name,ends_with("_P"))  %>%
  rename_at(.vars = vars(ends_with("_P")),
            .funs = funs(sub("_P$", "", .))) %>%
  column_to_rownames("name") %>%
  relocate(Serum,.after = "Human") %>%
  rename_with(~ ifelse(. == "Human", "Human\nserum", .)) %>%
  as.matrix()


FDR_df <- Plot_df %>%
  dplyr::select(name, ends_with("q")) %>%
  rename_at(.vars = vars(ends_with("_q")),
            .funs = funs(sub("_q$", "", .))) %>%
  column_to_rownames("name") %>%
  relocate(Serum,.after = "Human") %>%
  rename_with(~ ifelse(. == "Human", "Human\nserum", .)) %>%
  as.matrix()

FC_df[is.na(FC_df)] <- 0
P_df[is.na(P_df)] <- 1
FDR_df[is.na(FDR_df)] <- 1

p_1 <- P_df
p_1[P_df >= 0.05] <- 1
p_1[P_df < 0.05] <- 0
FDR_df <- FDR_df + p_1

class_level <- c("Lipid classes", "Amino acids and derivatives")
row_df1 <- data.frame(name = rownames(FC_df)) %>%
  mutate(Group = c(rep("Lipid classes",5),rep("Amino acids and derivatives",8))) %>%
  column_to_rownames("name") %>%
  dplyr::select(Group) %>%
  as.vector()
## Draw Heatmap ----------------------------------------------------------------
col_fun = circlize::colorRamp2(c(-2, 0, 2),  c("#01befe", "white", "#ff7d00"))
p <- Heatmap(FC_df, name = "Statistic", 
             col = col_fun, 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(P_df[i,j] <= 0.05) {
                 grid.points(x , y , pch = 1, size = unit(2, "mm"), gp = gpar(col = "black",lwd=2))
               }
               if(FDR_df[i,j] <= 0.25) {
                 grid.points(x , y , pch = 2, size = unit(4, "mm"), gp = gpar(col = "black",lwd=2))
               }
             },
             row_names_max_width = max_text_width(rownames(FC_df),gp = gpar(fontsize = 12)),
             column_names_max_height = max_text_width(colnames(FC_df),gp = gpar(fontsize = 12)),
             row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),
             row_split = factor(row_df1$Group,levels = class_level),
             cluster_rows = FALSE, cluster_columns = FALSE,
             column_names_side = "top",column_names_rot = 45, 
             show_row_names = T, #show_column_names = T,
             heatmap_legend_param = list(direction = "horizontal")
)
pdf("Mice_T_match_to_human.pdf",
    width = 6, height = 5)
draw(p, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()
