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
library(DoubletFinder)
library(ggplot2)
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

################################################################################ Loading data

# General function to process Seurat objects for any tissue
process_tissue_data <- function(data_dir, project_name, mt_pattern = "^mt-", min_cells = 3, min_features = 200, 
                                nFeature_min = 200, nFeature_max = 5000, percent_mt_max = 10, min_gene_count = 10) {
  tissue_data <- Read10X(data.dir = data_dir)
  tissue <- CreateSeuratObject(counts = tissue_data, project = project_name, min.cells = min_cells, min.features = min_features)
  tissue[["percent.mt"]] <- PercentageFeatureSet(tissue, pattern = mt_pattern)
  tissue <- subset(tissue, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < percent_mt_max)
  
  # Filter genes expressed in at least `min_gene_count` cells
  gene_counts <- rowSums(tissue@assays$RNA$counts > 0)
  selected_genes <- names(gene_counts[gene_counts >= min_gene_count])
  tissue <- tissue[selected_genes, ]

  tissue <- SCTransform(tissue, return.only.var.genes = FALSE)
  tissue <- RunPCA(tissue)
  tissue <- FindNeighbors(tissue, dims = 1:20)
  tissue <- FindClusters(tissue, resolution = 0.1)
  tissue <- RunUMAP(tissue, dims = 1:20)
  
  return(tissue)
}

# Liver
liver2 <- process_tissue_data(data_dir = "liver02/", project_name = "liver2")
liver6 <- process_tissue_data(data_dir = "liver06/", project_name = "liver6")
liver10 <- process_tissue_data(data_dir = "liver10/", project_name = "liver10")
liver14 <- process_tissue_data(data_dir = "liver14/", project_name = "liver14")
liver18 <- process_tissue_data(data_dir = "liver18/", project_name = "liver18")
liver22 <- process_tissue_data(data_dir = "liver22/", project_name = "liver22")


################################################################################ Doublet removal (using liver samples as examples)

# pK Identification (no ground-truth)
sweep.res.list_liver2 <- paramSweep(liver2, PCs = 1:20, sct = T)
sweep.stats_liver2 <- summarizeSweep(sweep.res.list_liver2, GT = FALSE)
bcmvn_liver2 <- find.pK(sweep.stats_liver2)
# Homotypic Doublet Proportion Estimate
annotations <- liver2@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver2@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# Run DoubletFinder with varying classification stringencies
liver2 <- doubletFinder(liver2, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver2$pK[which(bcmvn_liver2$BCmetric == max(bcmvn_liver2$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver2 <- subset(liver2, subset = DF.classifications_0.25_0.05_599 == "Singlet")

sweep.res.list_liver6 <- paramSweep(liver6, PCs = 1:20, sct = T)
sweep.stats_liver6 <- summarizeSweep(sweep.res.list_liver6, GT = FALSE)
bcmvn_liver6 <- find.pK(sweep.stats_liver6)
annotations <- liver6@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver6@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
liver6 <- doubletFinder(liver6, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver6$pK[which(bcmvn_liver6$BCmetric == max(bcmvn_liver6$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver6 <- subset(liver6, subset = DF.classifications_0.25_0.23_536 == "Singlet")


sweep.res.list_liver10 <- paramSweep(liver10, PCs = 1:20, sct = T)
sweep.stats_liver10 <- summarizeSweep(sweep.res.list_liver10, GT = FALSE)
bcmvn_liver10 <- find.pK(sweep.stats_liver10)
annotations <- liver10@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver10@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
liver10 <- doubletFinder(liver10, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver10$pK[which(bcmvn_liver10$BCmetric == max(bcmvn_liver10$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver10 <- subset(liver10, subset = DF.classifications_0.25_0.05_561 == "Singlet")

sweep.res.list_liver14 <- paramSweep(liver14, PCs = 1:20, sct = T)
sweep.stats_liver14 <- summarizeSweep(sweep.res.list_liver14, GT = FALSE)
bcmvn_liver14 <- find.pK(sweep.stats_liver14)
annotations <- liver14@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver14@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
liver14 <- doubletFinder(liver14, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver14$pK[which(bcmvn_liver14$BCmetric == max(bcmvn_liver14$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver14 <- subset(liver14, subset = DF.classifications_0.25_0.3_794 == "Singlet")

sweep.res.list_liver18 <- paramSweep(liver18, PCs = 1:20, sct = T)
sweep.stats_liver18 <- summarizeSweep(sweep.res.list_liver18, GT = FALSE)
bcmvn_liver18 <- find.pK(sweep.stats_liver18)
annotations <- liver18@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver18@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
liver18 <- doubletFinder(liver18, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver18$pK[which(bcmvn_liver18$BCmetric == max(bcmvn_liver18$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver18 <- subset(liver18, subset = DF.classifications_0.25_0.04_557 == "Singlet")

sweep.res.list_liver22 <- paramSweep(liver22, PCs = 1:20, sct = T)
sweep.stats_liver22 <- summarizeSweep(sweep.res.list_liver22, GT = FALSE)
bcmvn_liver22 <- find.pK(sweep.stats_liver22)
annotations <- liver22@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(liver22@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
liver22 <- doubletFinder(liver22, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn_liver22$pK[which(bcmvn_liver22$BCmetric == max(bcmvn_liver22$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
liver22 <- subset(liver22, subset = DF.classifications_0.25_0.01_723 == "Singlet")


################################################################################ Data merge and integration (using liver samples as examples)

DefaultAssay(liver2) <- "RNA"
DefaultAssay(liver6) <- "RNA"
DefaultAssay(liver10) <- "RNA"
DefaultAssay(liver14) <- "RNA"
DefaultAssay(liver18) <- "RNA"
DefaultAssay(liver22) <- "RNA"

liver <- merge(liver14, y = c(liver18, liver22, liver2, liver6, liver10), 
               add.cell.ids = c("liver14", "liver18", "liver22", "liver2", "liver6", "liver10"), project = "liver")

liver <- SCTransform(liver)
liver <- RunPCA(liver)
ElbowPlot(liver, ndims = 50)
liver <- FindNeighbors(liver, dims = 1:50)
liver <- FindClusters(liver, resolution = 0.1)
liver <- RunUMAP(liver, dims = 1:50)

liver <- IntegrateLayers(object = liver, method = CCAIntegration, normalization.method = "SCT", verbose = F)

liver <- FindNeighbors(liver, reduction = "integrated.dr", dims = 1:50)
liver <- FindClusters(liver, resolution = 0.8)
liver <- RunTSNE(liver, dims = 1:50, reduction = "integrated.dr")
liver <- RunUMAP(liver, dims = 1:50, reduction = "integrated.dr")

DimPlot(liver, group.by = c("orig.ident", "Phenotype", "SCT_snn_res.0.8"), reduction = "tsne")
DimPlot(liver, split.by = c("Phenotype"), reduction = "tsne", group.by = c("orig.ident", "SCT_snn_res.0.8"))


################################################################################ scType automatic cell type annotation (using liver samples as examples)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "liver_ref_combined.xlsx";
tissue = "Liver"

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = liver[["SCT"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(liver@meta.data$SCT_snn_res.0.8), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(liver@meta.data[liver@meta.data$SCT_snn_res.0.8==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(liver@meta.data$SCT_snn_res.0.8==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

liver@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  liver@meta.data$customclassif[liver@meta.data$SCT_snn_res.0.8 == j] = as.character(cl_type$type[1])
}
