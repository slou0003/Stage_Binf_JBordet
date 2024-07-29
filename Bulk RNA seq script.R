##############################################################################
# Title: RNA Sequencing Data Analysis using Seurat
# Author: Loutani Sara
# Date: 14/07/2024
# Institution: ULB Université Libre de Bruxelles _ MIU Lab Jules Bordet Institut
# 
# Description: 
# This script performs bulk RNA sequencing data analysis using the Seurat package.
# It includes steps for data loading, quality control, normalization, 
# dimensionality reduction, clustering, and visualization.
# 
# Requirements:
# - R version 4.0
# - Seurat package version [Seurat Version]
# - Additional required packages: dplyr, ggplot2, etc.
# 
# Usage:
# Adjust the file paths and parameters as needed.
# Run the script line by line or source it as a whole.
#
# Notes:
# - Ensure that the required packages are installed before running the script.
# - This script is intended for educational purposes and may require 
#   modifications for specific use cases.
##############################################################################





library(Seurat)
library(SeuratObject)
# Install openxlsx package if not already installed
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)


# Define file paths
counts_file <- "D:/Sara/Counts_MD_all.xlsx"
metadata_file <- "D:/Sara/metadataGCall_Cohort.xlsx"


### Read count table
library(openxlsx)
countData <- read.xlsx(counts_file)
rownames(countData) <- countData[,1]
countData <- countData[,-1]  # Remove the first column as it is now row names



# Read metadata
metadata <- read.xlsx(metadata_file)



### Create Seurat Object with Countdata
seurat= CreateSeuratObject(countData, project = "TLS_RNAseq", assay ="RNA")


#### Add Metadata to Seurat Object

dat1 = metadata$MD.status
dat2 = metadata$Cohort
dat3 = metadata$Samples
dat4 = metadata$Patients
dat5 = metadata$Whole.section.clusters
dat6 = metadata$TLS.Location
dat7 = metadata$TLS.surface.scored.µm2
  
seurat = AddMetaData(seurat, metadata = dat1, col.name = "MD.status")
seurat = AddMetaData(seurat, metadata = dat2, col.name = "Cohort")
seurat = AddMetaData(seurat, metadata = dat3, col.name = "Patients")
seurat = AddMetaData(seurat, metadata = dat4, col.name = "Samples")
seurat = AddMetaData(seurat, metadata = dat5, col.name = "Whole.section.clusters")
seurat = AddMetaData(seurat, metadata = dat6, col.name = "TLS.Location")
seurat = AddMetaData(seurat, metadata = dat7, col.name = "TLS.surface.scored.µm2")
seurat@meta.data




#### Normalization and finding variable features

seurat <- NormalizeData(seurat, normalization.method = "RC", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xlim = c(0, 30), ylim = c(0, 30))

print(plot1 + plot2)




#### Scaling data and running PCAScaling data and running PCA
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")

DimPlot(seurat, reduction = "pca")


DimHeatmap(seurat, dims = 1, cells = 66, balanced = TRUE)

DimHeatmap(seurat, dims = 1:15, cells = 66, balanced = TRUE)




#### Clustering and UMAP/TSNE
seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 1.2)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seurat <- RunUMAP(seurat, dims = 1:10)

DimPlot(seurat, reduction = "umap", dims = 1:2)  # Plotting dimensions 1 and 2
#to visualize high-dimensional data in a lower-dimensional space
seurat = RunTSNE(seurat, perplexity=15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

#The difference between the two DimPlot commands lies in the dimensionality reduction 
#techniques they use: t-SNE and UMAP. Both of these methods are commonly used in 
#single-cell RNA sequencing data analysis to visualize high-dimensional data in a 
#lower-dimensional space (typically 2D).

DimPlot(seurat, reduction = "tsne", label = T, label.size = 5)
DimPlot(seurat, reduction = "tsne", label = T, group.by = "Samples")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "MD.status")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "Cohort")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "Patients")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "Whole.section.clusters")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "TLS.Location")
DimPlot(seurat, reduction = "tsne", label = T, group.by = "TLS.surface.scored.µm2")
DimPlot(seurat, reduction = "umap", label = T, label.size = 5)

DimPlot(seurat, reduction = "umap", label = T, group.by = "Samples")
DimPlot(seurat, reduction = "umap", label = T, group.by = "MD.status")
DimPlot(seurat, reduction = "umap", label = T, group.by = "Cohort")
DimPlot(seurat, reduction = "umap", label = T, group.by = "Patients")
DimPlot(seurat, reduction = "umap", label = T, group.by = "Whole.section.clusters")
DimPlot(seurat, reduction = "umap", label = T, group.by = "TLS.Location")
DimPlot(seurat, reduction = "umap", label = T, group.by = "TLS.surface.scored.µm2")



# DimPlot(seurat, reduction = "umap", label = T, group.by = "Group")

DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "Samples" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "MD.status" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "Cohort" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "Patients" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "Whole.section.clusters" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "TLS.Location" )
DimPlot(seurat, reduction = "umap", label = T, label.size = 5, label.color = "black", split.by = "TLS.surface.scored.µm2" )




#### Saving and loading the Seurat object

saveRDS(seurat, file = "MD_seurat.rds")

seurat= readRDS(file = "MD_seurat.rds")






#### Finding markers and plotting

## All marekrs
library(dplyr)
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(seurat.markers, file = "seurat.markers.csv", sep=";", dec=",", row.names = T)

#to identify the top 50 markers for each cluster based on their average log fold change
top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(seurat, features = top10$gene, group.by = "MD.status")

DoHeatmap(seurat, features = markers, group.by = "MD.status")

#all_markers= seurat.markers[seurat.markers$p_val_adj < 0.05 & (is.finite(rowSums(seurat.markers))), ]
all_markers= seurat.markers[seurat.markers$p_val_adj < 0.05, ]

all_markers

write.table(all_markers, file = "all_markers-significant.csv", sep=";", dec=",", row.names = T)


all_markers_top50 <- all_markers %>% top_n(n = 50, wt = avg_log2FC)
all_markers_top50

markers.50= rownames(all_markers_top50)

DoHeatmap(seurat, features = markers.50)

DotPlot(seurat, features = markers.50, cols = c("blue", "red", "yellow", "green"), dot.scale = 8, group.by = "MD.status") +
  RotatedAxis()



#find conserved markers across conditions !!!!!!!!!!!!!!!!ERRROS
DefaultAssay(seurat) <- "RNA"
Idents(seurat) <- "RNA_snn_res.1.2"
table(seurat$RNA_snn_res.1.2, seurat$MD.status)

t= table(seurat$RNA_snn_res.1.2, seurat$MD.status, seurat$Samples, seurat$Cohort, seurat$Patients, seurat$TLS.Location, seurat$TLS.surface.scored.µm2)

write.table(t, file = "Repartition_samples_per_clusters.csv", sep= " ")


plot_heatmap(dataset = seurat, 
             markers = top10$gene,
             sort_var = c("seurat_clusters","Group"),
             anno_var = c("seurat_clusters","Group","Status"),
             anno_colors = list("Dark2",
                                c("red","orange","yellow","purple","blue","green"),
                                c("cadetblue", "aquamarine", "darkorchid1", "burlywood")),
             hm_limit = c(-1,0,1),
             hm_colors = c("purple","black","yellow"))

plot_heatmap(dataset = seurat, 
             markers = markers,
             sort_var = c("seurat_clusters","Group"),
             anno_var = c("seurat_clusters","Group","Status"),
             anno_colors = list("Dark2",
                                c("red","orange","yellow","purple","blue","green"),
                                c("cadetblue", "aquamarine", "darkorchid1", "burlywood")),
             hm_limit = c(-1,0,1),
             hm_colors = c("purple","black","yellow"))

# find all markers of cluster 0
cluster0.markers <- FindMarkers(seurat, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)

cluster0= cluster0.markers[cluster0.markers$p_val_adj < 0.05 & (is.finite(rowSums(cluster0.markers))), ]
cluster0
markers.to.plot = rownames(cluster0)
DoHeatmap(seurat, features = markers.to.plot, group.by = "MD.status") + NoLegend()
 
write.table(cluster0, file = "cluster0.markers_signif.csv", sep=";", dec=",", row.names = T)


cluster0_top50 <- cluster0 %>% top_n(n = 50, wt = avg_log2FC)
cluster0_top50 


markers.to.plot = rownames(cluster0_top50)
DotPlot(seurat, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()




# find all markers of cluster 1
cluster1.markers <- FindMarkers(seurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster1= cluster1.markers[cluster1.markers$p_val_adj < 0.05 & (is.finite(rowSums(cluster1.markers))), ]

write.table(cluster1, file = "cluster1.markers_signif.csv", sep=";", dec=",", row.names = T)


# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster2= cluster2.markers[cluster2.markers$p_val_adj < 0.05 & (is.finite(rowSums(cluster2.markers))), ]

write.table(cluster2, file = "cluster2.markers_signif.csv", sep=";", dec=",", row.names = T)



# find all markers of cluster 3 (fewer cells, so change number of ident)
cluster3.markers <- FindMarkers(seurat, ident.1 = 1, min.pct = 0.25)

head(cluster3.markers, n = 5)
cluster3= cluster3.markers[cluster3.markers$p_val_adj < 0.05 & (is.finite(rowSums(cluster3.markers))), ]

write.table(cluster3, file = "cluster3.markers_signif.csv", sep=";", dec=",", row.names = T)


# find all markers distinguishing one cluster from the others from 

### cluster0 = Mostly Active and Non-active

cluster0.vs_cluster1.2<- FindMarkers(seurat, ident.1 = 0, ident.2 = c(1,2), min.pct = 0.25)
head(cluster0.vs_cluster1.2, n = 5)
cluster0.vs_cluster1.2= cluster0.vs_cluster1.2[cluster0.vs_cluster1.2$p_val_adj < 0.05 & (is.finite(rowSums(cluster0.vs_cluster1.2))), ]


write.table(cluster0.vs_cluster1.2, file = "cluster0.vs_cluster1.2.markers.csv", sep=";", dec=",", row.names = T)




### cluster1 = Mostly Inactive aggreg and Non active

cluster1.vs_cluster0.2<- FindMarkers(seurat, ident.1 = 1, ident.2 = c(0,2), min.pct = 0.25)
head(cluster1.vs_cluster0.2, n = 5)
cluster1.vs_cluster0.2= cluster1.vs_cluster0.2[cluster1.vs_cluster0.2$p_val_adj < 0.05 & (is.finite(rowSums(cluster1.vs_cluster0.2))), ]

write.table(cluster1.vs_cluster0.2, file = "cluster1.vs_cluster0.2.markers.csv", sep=";", dec=",", row.names = T)





### cluster2 = A mix

cluster2.vs_cluster0.1<- FindMarkers(seurat, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
head(cluster2.vs_cluster0.1, n = 5)

cluster2.vs_cluster0.1= cluster2.vs_cluster0.1[cluster2.vs_cluster0.1$p_val_adj < 0.05 & (is.finite(rowSums(cluster2.vs_cluster0.1))), ]


write.table(cluster2.vs_cluster0.1, file = "cluster2.vs_cluster0.1.markers.csv", sep=";", dec=",", row.names = T)




### Cluster3 = Only non-active (already done in previous analysis)
# Compare cluster 3 vs clusters 0, 1, and 2
cluster3.vs_cluster0.1.2 <- FindMarkers(seurat, ident.1 = 3, ident.2 = c(0, 1, 2), min.pct = 0.25)
cluster3.vs_cluster0.1.2 <- cluster3.vs_cluster0.1.2[
  cluster3.vs_cluster0.1.2$p_val_adj < 0.05 & (is.finite(rowSums(cluster3.vs_cluster0.1.2))), ]
write.table(cluster3.vs_cluster0.1.2, file = "cluster3.vs_cluster0.1.2.markers.csv", sep = ";", dec = ",", row.names = TRUE)





# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# Load necessary libraries
library(Seurat)
library(dplyr)

# Assuming 'seurat' is your Seurat object

# Step 1: Find markers for each cluster compared to all remaining cells
# Load necessary libraries
library(Seurat)


# Step 1: Find markers for every cluster compared to all remaining cells
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Optionally, you can extract top markers for each cluster if needed
top_markers <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(seurat, features = c("CYP4F25P", "IGLV3-9", "IGKV3D-15", "IGLVI-70", "LALBA", "MUC19", "PRAMEF26", "PVRIG"))

# you can plot raw counts as well
VlnPlot(seurat, features = c("CYP4F25P", "IGLV3-9", "IGKV3D-15", "IGLVI-70", "LALBA", "MUC19", "PRAMEF26", "PVRIG"), slot = "counts", log = TRUE)

# Step 2: Find markers for each individual cluster
cluster0.markers <- FindMarkers(seurat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster1.markers <- FindMarkers(seurat, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster2.markers <- FindMarkers(seurat, ident.1 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster3.markers <- FindMarkers(seurat, ident.1 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Optionally, you can export the results to CSV files
write.csv(seurat.markers, file = "all_clusters_markers.csv", row.names = FALSE)
write.csv(cluster0.markers, file = "cluster0_markers.csv", row.names = FALSE)
write.csv(cluster1.markers, file = "cluster1_markers.csv", row.names = FALSE)
write.csv(cluster2.markers, file = "cluster2_markers.csv", row.names = FALSE)
write.csv(cluster3.markers, file = "cluster3_markers.csv", row.names = FALSE)



#### Saving and loading the Seurat object

saveRDS(seurat, file = "MD_seurat.rds")

seurat= readRDS(file = "MD_seurat.rds")

# heatmap top 8
top_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 8) %>%
  ungroup() -> top8
DoHeatmap(seurat, features = top8$gene) + NoLegend()


# Atribute cell types to clusters
new.cluster.ids <- c("B.0", "B.1", "B.2", "B.3")
# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(seurat)
pbmc_subset <- RenameIdents(breast_subset, selected_celltypes)
DimPlot(breast_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

saveRDS(seurat, file = "MD_seurat.rds")



############ STEP3 : DIFF GC+(A) et GC-(A ET B) see how they cluster differently
