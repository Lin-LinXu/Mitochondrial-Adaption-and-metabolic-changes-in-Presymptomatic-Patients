# Title: C13 Serine metabolite
# Project: Mitochondrial Adaption at the Center of Early Metabolic Changes in Presymptomatic Sepsis Patients
# Maintainer: Lin-Lin Xu
# Last update: Nov 28, 2024
# License: MIT

setwd("Codes/data")

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(magrittr)
library(ComplexHeatmap)
library(grid)


liver <- readRDS(file = "energy/liver_energy_subset.rds")
kidney <- readRDS(file = "energy/kidney_energy_subset.rds")
brain <- readRDS(file = "energy/brain_energy_subset.rds")
wat <- readRDS(file = "energy/wat_energy_subset.rds")

################################################################################ Liver hub genes ######################################################################
liver <- SetupForWGCNA(
  liver,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "liver_net" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
liver <- MetacellsByGroups(
  seurat_obj = liver,
  group.by = c("Cell_type", "SampleID"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on #integrated.dr
  k = 30, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Cell_type', # set the Idents of the metacell seurat object
  min_cells = 100,
  slot = "data",
  assay = 'SCT'
)

liver <- NormalizeMetacells(liver)

# set expression matrix for hdWGCNA
liver <- SetDatExpr(
  liver,
  group_name = c("Hepatocytes","Endothelial cells","Kupffer cells","Hepatic stellate cells",
                 "Cholangiocytes","B lymph cells", "Neutrophils"), # the name of the group of interest in the group.by column 
  group.by='Cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'SCT', # using SCT assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
liver <- TestSoftPowers(
  liver,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(liver)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(liver)
head(power_table)

# construct co-expression network:
liver <- ConstructNetwork(
  liver,
  soft_power = 10,
  overwrite_tom = TRUE,
  tom_name = 'liver_network' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(liver, main='liver cells energy hdWGCNA Dendrogram')
TOM <- GetTOM(liver)

# TOM was used as input to Cytoscape for network construction
# write.table(TOM, "liver_TOM.txt", sep = "\t")

liver <- ScaleData(liver, features=VariableFeatures(liver))
#liver <- ModuleEigengenes(liver)
liver <- ModuleEigengenes(
  liver,
  group.by.vars="SampleID"
)

# module eigengenes:
MEs <- GetMEs(liver, harmonized=T)

# compute eigengene-based connectivity (kME):
liver <- ModuleConnectivity(
  liver,
  group.by = 'Cell_type', 
  group_name = c("Hepatocytes","Endothelial cells","Kupffer cells","Hepatic stellate cells", "T lymph / NK cells",
                 "Cholangiocytes","B lymph cells", "Neutrophils"),
)

# rename the modules
liver <- ResetModuleNames(
  liver,
  new_name = "liver-M"
)

p <- PlotKMEs(liver, ncol=5)
p

# get the module assignment table:
modules <- GetModules(liver) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])

hub_liver <- GetHubGenes(liver, n_hubs = 10)
head(hub_liver)

################################################################################ Liver hub gene heatmap

all_cell_phe_energy_markers <- read.table("all_cell_phe_energy_markers.tsv", sep = "\t", header = T)
all_cell_phe_energy_markers <- all_cell_phe_energy_markers %>%
  left_join(hub_liver %>% select(gene_name, module), 
            by = c("Gene_symbol" = "gene_name"))

hub_liver <- all_cell_phe_energy_markers %>%
  filter(Gene_symbol %in% hub_liver$gene_name) %>%
  filter(tissue == "liver")

heatmap_data <- hub_liver %>%
  select(Gene_symbol, cell_type, avg_log2FC) %>%
  spread(key = cell_type, value = avg_log2FC) %>%
  as.data.frame()
heatmap_data[is.na(heatmap_data)] <- 0  
  
expanded_df <- heatmap_data %>%
  right_join(hub_liver %>% select(Gene_symbol = Gene_symbol, module), by = "Gene_symbol")

organs <- c("brain", "kidney", "liver", "wat")

heatmap_data <- expanded_df %>%
  select(-Gene_symbol, -module) %>%
  as.matrix()
rownames(heatmap_data) <- expanded_df$Gene_symbol

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
row_annotation <- rowAnnotation(Gene_symbol = anno_text(expanded_df$Gene_symbol, which = "row", gp = gpar(fontsize = 10)),
                                #Pathway = expanded_df$Pathway, 
                                category = expanded_df$module,
                                col = list(foo = col_fun))

row_split_gap = unit(1, "mm")
column_split_gap = unit(1, "mm")

row.names(expanded_df) <- expanded_df[,1]

heatmap_data <- data.frame(heatmap_data)

Heatmap(heatmap_data, name = "log2FC", 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 10), 
        row_dend_side = "left",
        right_annotation = row_annotation,
        cluster_rows = T, # 不进行聚类
        cluster_columns = T, # 不进行聚类
        #column_split = organ_info$Organ, # 按器官分组
        row_split = expanded_df$module,
        row_gap = row_split_gap, # 调整行分组之间的间距
        column_gap = column_split_gap, # 调整列分组之间的间距
        #row_split = expanded_df$Pathway, # 按代谢通路拆分
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        show_row_names = F
        #width = unit(100, "mm"),
        #height = unit(600, "mm")
)
  