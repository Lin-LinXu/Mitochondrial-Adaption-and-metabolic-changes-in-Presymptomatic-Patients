# Title: C13 Serine metabolite
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Zhengyuan Zhou
# Last update: Nov 28, 2024
# License: MIT

setwd("Codes/data")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(tibble)

################################################################################ DEGs at tissue level ###################################################################################

# Read energy metabolism-related genes
energy_genes <- read.table("energy_gene_mito_rm.txt", header = FALSE, stringsAsFactors = FALSE)$V1

liver <- readRDS("liver_SCT_CCA_meta.rds")
kidney <- readRDS("kidney_SCT_CCA_meta.rds")
brain <- readRDS("brain_SCT_CCA_meta.rds")
wat <- readRDS("wat_SCT_CCA_meta.rds")

ProcessTissuePhenotypeMarkers <- function(
    tissue, 
    ident_1 = "CLP", 
    group_by = "Phenotype", 
    logfc_threshold = 0.5, 
    p_val_adj_threshold = 0.05
) 
  {
  tissue <- PrepSCTFindMarkers(tissue)
  
  markers_all <- FindMarkers(
    tissue,
    min.cells.group = 1,
    min.cells.feature = 1,
    min.pct = 0,
    logfc.threshold = 0,
    only.pos = FALSE,
    ident.1 = ident_1,
    group.by = group_by
  )
  
  markers_filtered <- markers_all %>%
    mutate(
      abs_avg_log2FC = abs(avg_log2FC),
      direction = ifelse(avg_log2FC > 0, ident_1, "Sham")
    ) %>%
    filter(abs_avg_log2FC > logfc_threshold & p_val_adj < p_val_adj_threshold)
  
  return(markers_filtered)
}

liver_phe_markers <- ProcessTissuePhenotypeMarkers(liver)
kidney_phe_markers <- ProcessTissuePhenotypeMarkers(kidney)
brain_phe_markers <- ProcessTissuePhenotypeMarkers(brain)
wat_phe_markers <- ProcessTissuePhenotypeMarkers(wat)

all_tissue_phe_markers <- bind_rows(
  liver_phe_markers %>% mutate(tissue = "liver"),
  kidney_phe_markers %>% mutate(tissue = "kidney"),
  brain_phe_markers %>% mutate(tissue = "brain"),
  wat_phe_markers %>% mutate(tissue = "wat")
)

all_tissue_phe_markers <- all_tissue_phe_markers %>%
  mutate(Gene_symbol = sub("\\..*", "", rownames(all_tissue_phe_markers)))

all_genes <- unique(all_tissue_phe_markers$Gene_symbol)
common_genes <- intersect(all_genes, energy_genes)
all_tissue_phe_energy_markers <- all_tissue_phe_markers %>%
  filter(Gene_symbol %in% common_genes)

all_tissue_phe_energy_markers <- all_tissue_phe_energy_markers %>%
  rownames_to_column(var = "ID")

write.table(all_tissue_phe_energy_markers, "all_tissue_phe_energy_markers.tsv", sep = "\t", row.names = F, quote = F)


################################################################################ DEGs at cell type level ###############################################################################

ProcessCellPhenotypeMarkers <- function(tissue, ident_1 = "CLP", ident_2 = "Sham", group_by = "Phenotype", log2fc_threshold = 0.5, p_val_adj_threshold = 0.05) {
  cell_types <- unique(tissue@meta.data$Cell_type)
  deg_results <- list()

  for (cell_type in cell_types) {
    tissue_subset <- subset(tissue, subset = Cell_type == cell_type)
    tissue_subset <- PrepSCTFindMarkers(tissue_subset)
    
    markers <- FindMarkers(
      tissue_subset, 
      ident.1 = ident_1, 
      ident.2 = ident_2, 
      group.by = group_by,
      only.pos = FALSE, 
      recorrect_umi = FALSE
    )

    markers <- markers %>%
      mutate(abs_avg_log2FC = abs(avg_log2FC)) %>%
      filter(abs_avg_log2FC > log2fc_threshold & p_val_adj < p_val_adj_threshold) %>%
      mutate(
        direction = ifelse(avg_log2FC > 0, ident_1, ident_2),
        cell_type = cell_type
      ) %>%
      arrange(desc(abs_avg_log2FC))
    
    deg_results[[cell_type]] <- markers
  }
  
  all_deg_results <- bind_rows(deg_results)
  return(all_deg_results)
}

liver_deg <- ProcessCellPhenotypeMarkers(liver)
kidney_deg <- ProcessCellPhenotypeMarkers(kidney)
brain_deg <- ProcessCellPhenotypeMarkers(brain)
wat_deg <- ProcessCellPhenotypeMarkers(wat)

all_cell_phe_markers <- bind_rows(
  liver_deg %>% mutate(tissue = "liver"),
  kidney_deg %>% mutate(tissue = "kidney"),
  brain_deg %>% mutate(tissue = "brain"),
  wat_deg %>% mutate(tissue = "wat")
)

all_cell_phe_markers <- all_cell_phe_markers %>%
  mutate(Gene_symbol = sub("\\..*", "", rownames(all_cell_phe_markers)))

all_genes <- unique(all_cell_phe_markers$Gene_symbol)
common_genes <- intersect(all_genes, energy_genes)
all_cell_phe_energy_markers <- all_cell_phe_markers %>%
  filter(Gene_symbol %in% common_genes)

all_cell_phe_energy_markers <- all_cell_phe_energy_markers %>%
  rownames_to_column(var = "ID")

write.table(all_cell_phe_energy_markers, "all_cell_phe_energy_markers.tsv", sep = "\t", row.names = F, quote = F)
