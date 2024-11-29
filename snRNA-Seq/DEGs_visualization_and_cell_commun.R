# Title: C13 Serine metabolite
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

setwd("Codes/data")

library(ggplot2)
library(gridBase)
library(tidyr)
library(dplyr)
library(purrr)
library(tidyverse)
library(plotly)
library(ggrepel)
library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)

################################################################################ Energy-related gene extraction #######################################################
# Read energy metabolism-related genes from the text file
energy_genes <- read.table("energy_gene_mito_rm.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# Check the gene list
print(energy_genes)

# Function to extract relevant genes from a Seurat object
extract_genes_seurat <- function(seurat_obj, gene_list) {
  # Get all genes in the Seurat object (assume SCT assay is used)
  all_genes <- rownames(seurat_obj)
  
  # Find overlap between the gene list and genes in the Seurat object
  common_genes <- intersect(all_genes, gene_list)
  message("Number of genes found: ", length(common_genes))
  
  # Subset the Seurat object to keep only relevant genes
  seurat_subset <- subset(seurat_obj, features = common_genes)
  
  return(seurat_subset)
}

liver <- readRDS(file = "../data/liver_SCT_CCA_meta.rds")
kidney <- readRDS(file = "../data/kidney_SCT_CCA_meta.rds")
brain <- readRDS(file = "../data/brain_SCT_CCA_meta.rds")
wat <- readRDS(file = "../data/wat_SCT_CCA_meta.rds")

# Apply the function to each tissue-specific Seurat object
liver_subset <- extract_genes_seurat(liver, energy_genes)
kidney_subset <- extract_genes_seurat(kidney, energy_genes)
brain_subset <- extract_genes_seurat(brain, energy_genes)
wat_subset <- extract_genes_seurat(wat, energy_genes)

# Cell screen
FilterCellsByGeneExpression <- function(tissues, layer = "count", min_gene_counts = 1) {
  if (!is.list(tissues)) {
    stop("The input 'tissues' must be a list of single-cell objects.")
  }
  filtered_tissues <- list()
  for (tissue_name in names(tissues)) {
    tissue <- tissues[[tissue_name]]
    gene_expr <- GetAssayData(object = tissue, layer = layer)
    cell_gene_counts <- colSums(gene_expr > 0)
    filtered_cells <- which(cell_gene_counts >= min_gene_counts)
    tissue_filtered <- tissue[, filtered_cells]
    filtered_tissues[[tissue_name]] <- tissue_filtered
  }
  return(filtered_tissues)
}

tissues <- list(
  liver = liver_subset,
  kidney = kidney_subset,
  brain = brain_subset,
  wat = wat_subset
)

filtered_tissues <- FilterCellsByGeneExpression(tissues, layer = "count", min_gene_counts = 1)

liver <- filtered_tissues$liver
kidney <- filtered_tissues$kidney
brain <- filtered_tissues$brain
wat <- filtered_tissues$wat

# Save the resulting Seurat objects (optional)
saveRDS(liver, file = "../liver_energy_subset.rds")
saveRDS(kidney, file = "../kidney_energy_subset.rds")
saveRDS(brain, file = "../brain_energy_subset.rds")
saveRDS(wat, file = "../wat_energy_subset.rds")



################################################################################ Figure 5A ############################################################################
# Function to calculate average expression by phenotype
calculate_avg_expression <- function(seurat_obj, phenotype_column, tissue_name) {
  # Ensure the phenotype column exists
  if (!(phenotype_column %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Phenotype column", phenotype_column, "not found in metadata!"))
  }
  
  # Get SCT assay data (log-normalized counts assumed)
  expression_data <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  
  # Convert to a data frame for easier manipulation
  expression_df <- as.data.frame(t(expression_data))
  expression_df$Phenotype <- seurat_obj@meta.data[[phenotype_column]]
  
  # Calculate average expression for each phenotype and gene
  avg_expression_df <- expression_df %>%
    group_by(Phenotype) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(cols = -Phenotype, names_to = "Gene", values_to = "Expression") %>%
    mutate(Phenotype = paste(Phenotype, tissue_name, sep = "_")) %>%
    pivot_wider(names_from = Phenotype, values_from = Expression)
  
  return(avg_expression_df)
}

# Function to process all tissues and merge results
process_tissues <- function(tissue_objects, phenotype_column) {
  results <- list()
  
  for (tissue_name in names(tissue_objects)) {
    seurat_obj <- tissue_objects[[tissue_name]]
    message("Processing tissue: ", tissue_name)
    avg_expr <- calculate_avg_expression(seurat_obj, phenotype_column, tissue_name)
    results[[tissue_name]] <- avg_expr
  }
  
  # Merge all tissues into a single data frame by "Gene"
  combined_df <- Reduce(function(x, y) full_join(x, y, by = "Gene"), results)
  return(combined_df)
}

# Processing the four tissue-specific Seurat objects
tissues <- list(
  liver = liver,
  kidney = kidney,
  brain = brain,
  wat = wat
)

# Specify the phenotype column (e.g., "Phenotype")
phenotype_column <- "Phenotype"

# Calculate and combine average expression profiles
avg_expression_combined <- process_tissues(tissues, phenotype_column)

data_matrix <- avg_expression_combined

col_names <- colnames(data_matrix)
sham_columns <- grep("Sham", col_names, value = TRUE)
clp_columns <- grep("CLP", col_names, value = TRUE)
new_order <- c()
for (sham_col in sham_columns) {
  clp_col <- sub("Sham", "CLP", sham_col)
  new_order <- c(new_order, sham_col, clp_col)
}
data_matrix <- data_matrix[, new_order]

data_matrix_t <- t(data_matrix)
colnames(data_matrix_t) <- data_matrix_t[1,]
data_matrix_t <- data_matrix_t[-1,]
data_matrix_t[is.na(data_matrix_t)] <- 0
#str(data_matrix_t)
row <- rownames(data_matrix_t)
data_matrix_t <- as.data.frame(apply(data_matrix_t, 2, as.numeric))

pca_result <- prcomp(data_matrix_t, scale. = TRUE)
summary(pca_result)
pca_data <- as.data.frame(pca_result$x)
rownames(pca_data) <- row
pca_data$tissue <- rownames(pca_data)
head(pca_data)

ggplot(pca_data, aes(x = PC1, y = PC2, label = tissue, color = tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 1.5, vjust = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "PCA of Average Gene Expression Profiles of energy metabolism",
       x = "PC1",
       y = "PC2") +
  theme_minimal()



################################################################################ Figure 5B - Following the CellChat pipeline ##########################################
process_tissue_data <- function(tissue, data_path) {
  lapply(c("dplyr", "Seurat", "HGNChelper", "openxlsx", "CellChat"), library, character.only = TRUE)
  # Data input
  tissue_data <- readRDS(file.path(data_path, paste0(tissue, "_energy_subset.rds")))
  #tissue_meta <- read.table(file.path(meta_path, paste0(tissue, "_meta.tsv")), header = TRUE, row.names = 1, sep = "\t")
  #tissue_data <- AddMetaData(tissue_data, metadata = tissue_meta)
  
  # Energy gene screening
  #gene_list <- readLines(gene_list_path)
  #tissue_data <- tissue_data[gene_list, ]
  
  # Cell filtering
  gene_expr <- GetAssayData(object = tissue_data, layer = "count")
  cell_gene_counts <- colSums(gene_expr > 0)
  filtered_cells <- which(cell_gene_counts > 0)
  tissue_data <- tissue_data[, filtered_cells]
  
  # Grouped by phenotype
  tissue_CLP <- subset(tissue_data, subset = Phenotype == "CLP")
  tissue_sham <- subset(tissue_data, subset = Phenotype == "Sham")
  
  # CellChat
  cellchat_CLP <- createCellChat(object = tissue_CLP, group.by = "Cell_type", assay = "SCT")
  cellchat_sham <- createCellChat(object = tissue_sham, group.by = "Cell_type", assay = "SCT")
  
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- CellChatDB
  
  cellchat_CLP@DB <- CellChatDB.use
  cellchat_sham@DB <- CellChatDB.use
  
  cellchat_CLP <- subsetData(cellchat_CLP)
  cellchat_sham <- subsetData(cellchat_sham)
  
  future::plan("multisession", workers = 4)
  cellchat_CLP <- identifyOverExpressedGenes(cellchat_CLP)
  cellchat_CLP <- identifyOverExpressedInteractions(cellchat_CLP)
  cellchat_sham <- identifyOverExpressedGenes(cellchat_sham)
  cellchat_sham <- identifyOverExpressedInteractions(cellchat_sham)
  
  cellchat_CLP <- computeCommunProb(cellchat_CLP, type = "truncatedMean", trim = 0.05, raw.use = TRUE, population.size = FALSE)
  cellchat_sham <- computeCommunProb(cellchat_sham, type = "truncatedMean", trim = 0.05, raw.use = TRUE, population.size = FALSE)
  
  cellchat_CLP <- computeCommunProbPathway(cellchat_CLP)
  cellchat_sham <- computeCommunProbPathway(cellchat_sham)
  
  cellchat_CLP <- aggregateNet(cellchat_CLP)
  cellchat_sham <- aggregateNet(cellchat_sham)
  
  cellchat_CLP <- netAnalysis_computeCentrality(cellchat_CLP, slot.name = "netP")
  cellchat_sham <- netAnalysis_computeCentrality(cellchat_sham, slot.name = "netP")
  
  object.list <- list(sham = cellchat_sham, CLP = cellchat_CLP)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
  
  return(cellchat)
}

data_path <- "energy/"
#meta_path <- "liver"
#gene_list_path <- "Figures/new_energy_metabolism/mito_rm/energy_gene_mito_rm.txt"

cellchat_liver <- process_tissue_data("liver", data_path)
cellchat_kidney <- process_tissue_data("kidney", data_path)
cellchat_brain <- process_tissue_data("brain", data_path)
cellchat_wat <- process_tissue_data("wat", data_path)




################################################################################ Figure 5C - Common DEGs #######################################################
all_tissue_phe_energy_markers <- read.table("all_tissue_phe_energy_markers.tsv", sep = "\t", header = T)

all_tissue_phe_energy_markers$tissue <- factor(all_tissue_phe_energy_markers$tissue, levels = c("liver", "kidney", "brain", "wat"))
all_tissue_phe_energy_markers$direction <- factor(all_tissue_phe_energy_markers$direction, levels = c("Sham", "CLP"))

# Count the number of DEGs at the tissue level
summary_data <- all_tissue_phe_energy_markers %>%
  group_by(tissue, direction) %>%
  summarise(Gene_count = n(), .groups = "drop")

# Bar plot
ggplot(summary_data, aes(x = direction, y = Gene_count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(~ tissue, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("CLP" = "#e7211a", "Sham" = "#194591")) +
  labs(
    title = "DEGs at the tissue level",
    x = "Phenotype",
    y = "Number of DEGs"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


################################################################################ Figure 5D - ##########################################################################
# Identify DEGs with consistent direction across all tissues
common_DEGs <- all_tissue_phe_energy_markers %>%
  group_by(Gene_symbol) %>%
  summarise(tissue_count = n_distinct(tissue), direction_count = n_distinct(direction), .groups = "drop") %>%
  filter(tissue_count == 4 & direction_count == 1) %>%
  left_join(all_tissue_phe_energy_markers, by = "Gene_symbol")

all_cell_phe_energy_markers <- read.table("all_cell_phe_energy_markers.tsv", sep = "\t", header = T)
common_genes <- unique(common_DEGs$Gene_symbol)
t_DEGs_energy_markers <- all_cell_phe_energy_markers %>%
  filter(Gene_symbol %in% common_genes)

t_DEGs_energy_markers <- t_DEGs_energy_markers %>%
  mutate(cell_name = paste(cell_type, tissue, sep = "_"))

# 将数据转换为适合绘制热图的格式
heatmap_data <- t_DEGs_energy_markers %>%
  select(Gene_symbol, cell_name, avg_log2FC) %>%
  spread(key = cell_name, value = avg_log2FC) %>%
  as.data.frame()
heatmap_data[is.na(heatmap_data)] <- 0
#row.names(heatmap_data) <- heatmap_data[,1]
#heatmap_data <- heatmap_data[,-1]
#write.table(heatmap_data, "cell_level/common_energy_DEGs_cell_matrix.txt", sep = "\t")

pathway_info <- data.frame(read.table("common_DEGs/common_energy_pathway_info.txt", sep = "\t", header = T))
organ_info <- data.frame(read.table("common_DEGs/organ_info.txt", sep = "\t", header = T))
row.names(organ_info) <- organ_info[,1]
organ_info <- organ_info[,-1, drop = F]

#merged_df <- read.table("common_energy_DEGs_sham_matrix.txt", sep = "\t", header = T)

# 将基因-代谢通路数据处理为长格式
gene_pathways_long <- pathway_info %>%
  separate_rows(Pathway, sep = ";") %>%
  unique()

# 扩展表达矩阵，处理基因重复出现的情况
expanded_df <- heatmap_data %>%
  right_join(gene_pathways_long, by = "Gene_symbol")

# 对不同器官的数据分别进行标准化
organs <- c("brain", "kidney", "liver", "wat")

# 准备热图数据
heatmap_data <- expanded_df %>%
  select(-Gene_symbol, -Pathway) %>%
  as.matrix()
rownames(heatmap_data) <- expanded_df$Gene_symbol

# 准备行注释
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
row_annotation <- rowAnnotation(Gene_symbol = anno_text(expanded_df$Gene_symbol, which = "row", gp = gpar(fontsize = 10)),
                                #Pathway = expanded_df$Pathway, 
                                category = expanded_df$Pathway,
                                col = list(foo = col_fun))

# 准备分组信息
organ_groups <- rep(organs, each = 1) # each = 1

row_split_gap = unit(3, "mm")  # 设置行分组间的间距为4毫米

# 绘制热图并按代谢通路名称拆分
column_split_gap = unit(3, "mm")  # 设置分组间的间距为4毫米

row.names(pathway_info) <- pathway_info[,1]
pathway_info <- data.frame(pathway_info[,-1, drop = FALSE])

Heatmap(heatmap_data, name = "log2FC", 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 10), 
        row_dend_side = "left",
        right_annotation = row_annotation,
        cluster_rows = T, # 不进行聚类
        cluster_columns = T, # 不进行聚类
        column_split = organ_info$Organ, # 按器官分组
        #row_gap = row_split_gap, # 调整行分组之间的间距
        column_gap = column_split_gap, # 调整列分组之间的间距
        #row_split = expanded_df$Pathway, # 按代谢通路拆分
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        show_row_names = F,
        #width = unit(10, "mm"),
        #height = unit(6, "mm")
)
