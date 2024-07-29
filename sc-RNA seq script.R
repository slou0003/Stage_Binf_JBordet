##############################################################################
# Title: RNA Sequencing Data Analysis using Seurat
# Author: Loutani Sara
# Date: 14/07/2024
# Institution: ULB Université Libre de Bruxelles _ MIU Lab Jules Bordet Institut
# 
# Description: 
# This script performs RNA sequencing data analysis using the Seurat package.
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
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# Article's Processed datset path
Data <- "C:/Users/moisara/Desktop/Internship_Bordet/Blue_print_train/Data/panB_scRNA_processed_data.rds"

# Seurat objcet
data.obj <- readRDS(Data)


########### Split the seurat onject BY cancer TYPE :
########### I. Breast cancer BC subset :

###### Split the Seurat object by site type 
b_cell_subsets <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Breast", ]))

#head(b_cell_subsets@meta.data, 5)
Idents(b_cell_subsets) = "celltype"

#1.Vlnplot:
VlnPlot(b_cell_subsets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Features are unique genes 
# Ncounts are the expression amount
# percent is metocondrial <<  indiquates how alive the sample is

# 2.FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(b_cell_subsets, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(b_cell_subsets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Idents(data.obj) = "celltype"
# UMAPPlot(data.obj)

# Plot PCA
DimPlot(b_cell_subsets, reduction = "pca") + NoLegend()

# Heatmap
DimHeatmap(b_cell_subsets, dims = 1, cells = 500, balanced = TRUE)


# Breast cancer UMAP plot
UMAPPlot(b_cell_subsets)
#head(Idents(b_cell_subsets), 5)

# Dotplot all clusters all genes 
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Define your markers
B_marker <- c("TCL1A", "FCER2", "IL4R", "IGHD", #### NaiveB
              "IFIT3", "IFI44L", "STAT1", "ISG15", ### IFN
              "HSPA1A", "DNAJB1", ### Activated
              "MT1X", "MT2A", "SLC30A1",
              "EGR1", "DUSP2", #### ACB1
              "NR4A2", "CD69", "CD83", #### ACB2
              "CCR7", "PIM3", "NFKBID",
              "S100A10", "CRIP1", "S100A4", "ITGB1", "CD27", "CR2", "AIM2", "GPR183", "CD1C",
              "DUSP4", "FCRL5", "ZEB2", "ITGAX", "FGR", "FCRL4", "CD274",
              "NME1", "PSME2", "ENO1", "FABP5", ### PreGC
              "CXCR4", "GNB2L1", "ATP5L", "SUGCT", ### DZGC and GC
              "LMO2", "GMDS", "PRPSAP2", "MARCKSL1", ### LZGC
              "STMN1", "TUBB", "HMGB2", "TUBA1B", "MKI67","SARAF","CD79A","RPL18A","RPS27A","CD74","RPS19",
              "JCHAIN", "PRDM1", "XBP1", "MZB1")

hypoxia <- as.character(B_marker)

# Create Dot Plot
dot_plot <- DotPlot(object = b_cell_subsets, features = B_marker, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = hypoxia, labels = hypoxia)

# Display the Dot Plot
print(dot_plot)


# Extract and print all cluster names
cluster_names <- levels(Idents(b_cell_subsets))
print(cluster_names)



####### find all markers of cluster 1:
cluster1.markers <- FindMarkers(b_cell_subsets, ident.1 = "B.01.TCL1A+naiveB", min.pct = 0.25)
head(cluster1.markers, n = 20)

#Top 10 Expressed Genes in Cluster B.01.TCL1A+naiveB
BC1_top10_genes <- head(row.names(cluster1.markers), 10)

# Plot the top 10 expressed gene
VlnPlot(b_cell_subsets, features = BC1_top10_genes, pt.size = 1.5) 
# Print the plot
print(plot)

# plot the 2 first genes in Vlnplot
VlnPlot(b_cell_subsets, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2], row.names(cluster1.markers)[4]))

# find all markers distinguishing cluster 15 from clusters 1 and 14
cluster15.markers <- FindMarkers(b_cell_subsets, ident.1 = "B.15.Plasma cell", ident.2 = c("B.01.TCL1A+naiveB", "B.14.Plasmablast"), min.pct = 0.25)
head(cluster15.markers, n = 5)
VlnPlot(b_cell_subsets, features = c(row.names(cluster15.markers)[1], row.names(cluster15.markers)[2]))


####### ALL CLUSTERS COMPARAISON :

# find 10 top markers for every cluster compared to all remaining cells, report only the positive ones
## Find all markers 
b_subsets.markers <- FindAllMarkers(b_cell_subsets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10_genes_clusters <- head(row.names(b_subsets.markers), 10)
x <- top10_genes_clusters %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
# top << p value = 0
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(b_cell_subsets, features = x$gene[1:4])
FeaturePlot(b_cell_subsets, features = x$gene[5:8])
FeaturePlot(b_cell_subsets, features = x$gene[9:10])


# Do this for a given genes as a combined plot
p <- FeaturePlot(b_cell_subsets, features = top10_genes_clusters, combine = TRUE)
p <- lapply(X = p, FUN = function(x) x + 
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5)) +
              theme(axis.title.x = element_text(size = 5)) +
              theme(axis.text.y = element_text(size = 5)) +
              theme(axis.text.x = element_text(size = 5)) +
              theme(legend.position = "none")  )

CombinePlots(plots = p)

# Essaie find markers for all clusters , extract top 6 then plot:
# Find 6 markers for every cluster compared to all remaining cells, report only the positive ones

# Get the top 6 marker genes for all clusters
top6_genes_clusters <- head(row.names(b_subsets.markers), 6)

# FeaturePlot for each set of genes
p1 <- FeaturePlot(b_cell_subsets, features = top6_genes_clusters[1:4], combine = FALSE)
p2 <- FeaturePlot(b_cell_subsets, features = top6_genes_clusters[5:6], combine = FALSE)

# Adjust the themes for each plot
p1 <- lapply(p1, function(x) x + theme(plot.title = element_text(size = 8),
                                       axis.title.y = element_text(size = 5),
                                       axis.title.x = element_text(size = 5),
                                       axis.text.y = element_text(size = 5),
                                       axis.text.x = element_text(size = 5),
                                       legend.position = "none"))
p2 <- lapply(p2, function(x) x + theme(plot.title = element_text(size = 8),
                                       axis.title.y = element_text(size = 5),
                                       axis.title.x = element_text(size = 5),
                                       axis.text.y = element_text(size = 5),
                                       axis.text.x = element_text(size = 5),
                                       legend.position = "none"))


# Combine the plots
combined_plot <- CombinePlots(plots = list(p1, p2))

# Print the combined plot
print(combined_plot)


### TOP9 in each cluster
top9 <- b_subsets.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 9, wt = avg_log2FC)

top9


# noo Heatmap of top 9
p2 <- DoHeatmap(b_cell_subsets, features = top9$gene, group.bar.height = 0.01,size=3,combine = FALSE) 

p2 <- lapply(X = p2, FUN = function(x) x + 
               theme(plot.title = element_text(size = 8)) +
               theme(axis.title.y = element_text(size = 5)) +
               theme(axis.title.x = element_text(size = 5)) +
               theme(axis.text.y = element_text(size = 3)) +
               theme(legend.position = "none")  )

CombinePlots(plots = p2)


          
# Atribute cell types to clusters
new.cluster.ids <- c("B.01.TCL1A+naiveB", "B.03.HSP+B", "B.08.ITGB1+SwBm", "B.06.NR4A2+ACB2", "B.05.EGR1+ACB", 
                     "B.10.ENO1+Pre_GCB", "B.02.IFIT3+B", "B.07.CCR7+ACB3", "B.09.DUSP4+AtM", "B.14.Plasmablast", 
                     "B.15.Plasma cell", "B.12.LMO2+LZ_GCB", "B.04.MT1X+B", "B.11.SUGCT+DZ_GCB", "B.13.Cycling_GCB")
names(new.cluster.ids) <- levels(b_cell_subsets)

b_cell_subsets <- RenameIdents(b_cell_subsets, new.cluster.ids)
DimPlot(b_cell_subsets, reduction = "pca", label = TRUE, pt.size = 0.5)
DimPlot(b_cell_subsets, reduction = "umap", label = TRUE, pt.size = 0.5)

# Dotplot for all clusters of top 6 expressed genes in all clusters 

B_marker= c("TCL1A","FCER2","IL4R","IGHD" , ####NaiveB
            "IFIT3","IFI44L","STAT1","ISG15",###IFN
            "HSPA1A","DNAJB1",###Activated
            "MT1X","MT2A","SLC30A1",
            "EGR1","DUSP2",####ACB1
            "NR4A2","CD69","CD83",####ACB2
            "CCR7","PIM3","NFKBID",
            "S100A10","CRIP1","S100A4","ITGB1","CD27","CR2","AIM2","GPR183","CD1C",
            "DUSP4","FCRL5","ZEB2","ITGAX","FGR","FCRL4","CD274",
            "NME1","PSME2","ENO1","FABP5",###PreGC
            "CXCR4","GNB2L1","ATP5L","SUGCT",###DZGC and GC
            "LMO2","GMDS","PRPSAP2","MARCKSL1",###LZGC
            "STMN1","TUBB","HMGB2","TUBA1B","MKI67",
            "JCHAIN","PRDM1","XBP1","MZB1"
            
)
# Define the hypoxia variable
hypoxia <- top6_genes_clusters # Assuming hypoxia labels correspond to your genes

DotPlot(object = b_cell_subsets, features =  top6_genes_clusters,scale = T,group.by = "celltype") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks= hypoxia,labels=  hypoxia)

# Dotplot for all clusters of top 10 expressed genes

B_marker= c("TCL1A","FCER2","IL4R","IGHD" , ####NaiveB
            "IFIT3","IFI44L","STAT1","ISG15",###IFN
            "HSPA1A","DNAJB1",###Activated
            "MT1X","MT2A","SLC30A1",
            "EGR1","DUSP2",####ACB1
            "NR4A2","CD69","CD83",####ACB2
            "CCR7","PIM3","NFKBID",
            "S100A10","CRIP1","S100A4","ITGB1","CD27","CR2","AIM2","GPR183","CD1C",
            "DUSP4","FCRL5","ZEB2","ITGAX","FGR","FCRL4","CD274",
            "NME1","PSME2","ENO1","FABP5",###PreGC
            "CXCR4","GNB2L1","ATP5L","SUGCT",###DZGC and GC
            "LMO2","GMDS","PRPSAP2","MARCKSL1",###LZGC
            "STMN1","TUBB","HMGB2","TUBA1B","MKI67",
            "JCHAIN","PRDM1","XBP1","MZB1"
            
)
# Define the hypoxia variable
hypoxia <- top10_genes_clusters # Assuming hypoxia labels correspond to your genes

DotPlot(object = b_cell_subsets, features =  top10_genes_clusters,scale = T,group.by = "celltype") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks= hypoxia,labels=  hypoxia)


# Step2: Define the top 10 expressed genes for each cluster
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(b_cell_subsets, features = top10_genes_clusters, ncol = 3)
RidgePlot(b_cell_subsets, features = top10_genes_clusters, ncol = 3, combine = TRUE) +
  theme(axis.text.y = element_text(size = 3))

# Assuming top6_genes_clusters is defined and contains the top 6 marker genes for each cluster

# Load the ggridges package
library(ggridges)

# Create ridge plots for top 6 marker genes in each cluster
RidgePlot(b_cell_subsets, features = top6_genes_clusters, ncol = 2, combine = TRUE) +
  theme(axis.text.y = element_text(size = 6))


#### Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(b_cell_subsets, features = top6_genes_clusters)

##### check top 10 expressed genes in each cluster
# Load necessary libraries
library(dplyr)
library(Seurat)

####### Extract the top 10 markers of all the clusters 
###### Define the cluster names of interest
clusters_of_interest <- c("B.01.TCL1A+naiveB", "B.03.HSP+B", "B.08.ITGB1+SwBm", "B.06.NR4A2+ACB2", "B.05.EGR1+ACB",
                          "B.10.ENO1+Pre_GCB", "B.02.IFIT3+B", "B.07.CCR7+ACB3", "B.09.DUSP4+AtM", "B.14.Plasmablast",
                          "B.15.Plasma cell", "B.12.LMO2+LZ_GCB", "B.04.MT1X+B", "B.11.SUGCT+DZ_GCB", "B.13.Cycling_GCB")

# Function to get top 10 markers for a given cluster
get_top_markers <- function(cluster_name) {
  top_markers <- b_subsets.markers %>%
    filter(cluster == cluster_name) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(desc(avg_log2FC)) %>%
    pull(gene)
  return(top_markers)
}

# Loop through each cluster and print the top 10 markers
for (cluster_name in clusters_of_interest) {
  top_markers <- get_top_markers(cluster_name)
  cat("Top 10 markers for cluster", cluster_name, ":\n")
  print(top_markers)
  cat("\n") # Print a new line for better readability
}

##### TO know the cancers sites in seurat obj
# Extract unique cancer sites from metadata
cancer_sites <- unique(data.obj@meta.data$site)

# Print the unique cancer sites
print(cancer_sites)

################ II. Thyroid
###### Split the Seurat object by site type : Thyroide
thr_obj <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Thyroid", ]))
VlnPlot(thr_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(thr_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(thr_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

########### Umap THRC
Idents(thr_obj) = "celltype"
UMAPPlot(thr_obj)

# Plot PCA
DimPlot(thr_obj, reduction = "pca") + NoLegend()

# Heatmap
DimHeatmap(thr_obj, dims = 1, cells = 500, balanced = TRUE)
# THRc UMAP plot
UMAPPlot(thr_obj)

#DOTPLOT for all markers
# Define your markers
B_marker_thr <- c("TCL1A", "FCER2", "IL4R", "IGHD", #### NaiveB
              "IFIT3", "IFI44L", "STAT1", "ISG15", ### IFN
              "HSPA1A", "DNAJB1", ### Activated
              "MT1X", "MT2A", "SLC30A1",
              "EGR1", "DUSP2", #### ACB1
              "NR4A2", "CD69", "CD83", #### ACB2
              "CCR7", "PIM3", "NFKBID",
              "S100A10", "CRIP1", "S100A4", "ITGB1", "CD27", "CR2", "AIM2", "GPR183", "CD1C",
              "DUSP4", "FCRL5", "ZEB2", "ITGAX", "FGR", "FCRL4", "CD274",
              "NME1", "PSME2", "ENO1", "FABP5", ### PreGC
              "CXCR4", "GNB2L1", "ATP5L", "SUGCT", ### DZGC and GC
              "LMO2", "GMDS", "PRPSAP2", "MARCKSL1", ### LZGC
              "STMN1", "TUBB", "HMGB2", "TUBA1B", "MKI67","SARAF","CD79A","RPL18A","RPS27A","CD74","RPS19",
              "JCHAIN", "PRDM1", "XBP1", "MZB1")

hypoxia <- as.character(B_marker_thr)

# Create Dot Plot
dot_plot <- DotPlot(object = thr_obj, features = B_marker_thr, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = hypoxia, labels = hypoxia)

# Display the Dot Plot
print(dot_plot)

############ Markers comparaisons: TOP 10 in all clusters

# Find markers 
thr_subsets.markers <- FindAllMarkers(thr_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Find 10 markers for every cluster compared to all remaining cells, report only the positive ones
# Get the top 10 marker genes for all clusters
top10_genes_thr <- head(rownames(thr_subsets.markers), 10)

# FeaturePlot for each set of genes
p1 <- FeaturePlot(thr_obj, features = top10_genes_thr[1:10], pt.size = 0.5, min.cutoff = 0, combine = FALSE)


# Adjust the themes for each plot
p1 <- lapply(p1, function(x) x + theme(plot.title = element_text(size = 8),
                                       axis.title.y = element_text(size = 5),
                                       axis.title.x = element_text(size = 5),
                                       axis.text.y = element_text(size = 5),
                                       axis.text.x = element_text(size = 5),
                                       legend.position = "none"))

library(patchwork)
# Ensure p1 and p2 are not in list format
p1 <- p1[[1]]  # Extract the first plot from the list
p2 <- p2[[1]]  # Extract the second plot from the list

# Combine the plots using the patchwork system
combined_plot <- p1 + p2

# Print the combined plot
print(combined_plot)
















##### ETAPE 2 / 20-06-24: SUBSET SEURAT OBJECT BY SITE and by clsuster crteria B01 B08 AND B09

# Load necessary library
library(Seurat)

# Define the cell types to include
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")

# Colon cancer subset
colon_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Colon" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(colon_subset) <- "celltype"


# Liver cancer subset
liver_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Liver" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(liver_subset) <- "celltype"

# Gallbladder cancer subset
gallbladder_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Gallbladder" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(gallbladder_subset) <- "celltype"

# NB (Neuroblastoma) cancer subset
nb_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "NB" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(nb_subset) <- "celltype"

# PBMC (Peripheral Blood Mononuclear Cells) cancer subset
pbmc_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "PBMC" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(pbmc_subset) <- "celltype"

# Breast cancer subset
breast_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Breast" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(breast_subset) <- "celltype"

# Esophagus  subset
esophagus_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Esophagus" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(esophagus_subset) <- "celltype"

# HNSC (Head and Neck Squamous Cell Carcinoma) subset
hnsc_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "HNSC" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(hnsc_subset) <- "celltype"

# Kidney cancer subset
kidney_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Kidney" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(kidney_subset) <- "celltype"

# Skin cancer subset
skin_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Skin" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(skin_subset) <- "celltype"

# Lung cancer subset
lung_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Lung" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(lung_subset) <- "celltype"

# Endometrioid cancer subset
endometrioid_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Endometrioid" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(endometrioid_subset) <- "celltype"

# Pancreas cancer subset
pancreas_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Pancreas" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(pancreas_subset) <- "celltype"

# Thyroid cancer subset
thyroid_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Thyroid" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(thyroid_subset) <- "celltype"

# Bladder cancer subset
bladder_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Bladder" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(bladder_subset) <- "celltype"

# Cervix cancer subset
cervix_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Cervix" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(cervix_subset) <- "celltype"

# Ovary cancer subset
ovary_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Ovary" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(ovary_subset) <- "celltype"

# Stomach cancer subset
stomach_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "Stomach" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(stomach_subset) <- "celltype"

# GIST (Gastrointestinal Stromal Tumor) cancer subset
gist_subset <- subset(data.obj, cells = rownames(data.obj@meta.data[data.obj@meta.data$site == "GIST" & data.obj@meta.data$celltype %in% selected_celltypes, ]))
Idents(gist_subset) <- "celltype"





############################ 1- Breat cancer



# 1- Plot PCA
DimPlot(breast_subset, reduction = "pca") + NoLegend()

# 2- PCA
UMAPPlot(breast_subset)


# 3- Heatmap
DimHeatmap(breast_subset, dims = 1, cells = 500, balanced = TRUE)
   
# 5- Dotplot
# Define markers of intrest
B_markers <- c("AICDA", "CD79A", "CD79B", "FCRL5", "CD19", "MS4A1", "CD38", "MZB1", "IKZF3", "FCRL4", "FCRL3", "FCRLA")

# Create Dot Plot
dot_plot <- DotPlot(object = breast_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)



        
# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(breast_subset)
pbmc_subset <- RenameIdents(breast_subset, selected_celltypes)
DimPlot(breast_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

# Create Ridge Plot for each marker
RidgePlot(breast_subset, features = B_markers, ncol = 2) 

library(cowplot)

ridge_plots <- RidgePlot(breast_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)


# Violin plot - Visualize single cell expression distributions in each cluster : All cells have the same value of TCF21.
VlnPlot(breast_subset, features = B_markers)


# Feature plot - visualize feature expression in low-dimensional space: All cells have the same value (0) of “TCF21”
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(breast_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Breast cancer Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)


saveRDS(breast_subset, file = "BC_subset.rds")
## Fichier sur flash-disque




################################ 2- COLON CANCER


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)


# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(colon_subset)
DimPlot(colon_subset, reduction = "umap", label = TRUE, pt.size = 0.5)


# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(colon_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(colon_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(colon_subset, dims = 1, cells = 500, balanced = TRUE)

# 4. Create Dot Plot
dot_plot <- DotPlot(object = colon_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
ridge_plots <- RidgePlot(colon_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster
VlnPlot(colon_subset, features = B_markers)

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space
FeaturePlot(colon_subset, features = B_markers)

# Optional: Adjust the themes for each plot (if necessary)
# FeaturePlots for each marker 
library(patchwork)
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.5  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(colon_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)



###################### 3- Lymphe Node cancer

# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(ln_subset)

ln_subset <- RenameIdents(ln_subset, selected_celltypes)
DimPlot(ln_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(ln_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(ln_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(ln_subset, dims = 1, cells = 500, balanced = TRUE)

# 4. Create Dot Plot
dot_plot <- DotPlot(object = ln_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
ridge_plots <- RidgePlot(ln_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster:All cells have the same value of TCF21,All cells have the same value of CD8B2.
VlnPlot(ln_subset, features = B_markers)

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space: 1: All cells have the same value (0) of “CD8B2” 
#2: All cells have the same value (0) of “TCF21” 
FeaturePlot(ln_subset, features = B_markers)

# Optional: Adjust the themes for each plot (if necessary)
# FeaturePlots for each marker (optional)
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(ln_subset, features = marker, pt.size = 0.5, min.cutoff = 0, combine = FALSE) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "none")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3)
print(combined_plots)







############# 4- Liver cancer


# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(liver_subset)

liver_subset <- RenameIdents(liver_subset, selected_celltypes)
DimPlot(liver_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(liver_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(liver_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(liver_subset, dims = 1, cells = 500, balanced = TRUE)

# 4. Create Dot Plot
dot_plot <- DotPlot(object = liver_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
ridge_plots <- RidgePlot(liver_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster: All cells have the same value of TCF21
VlnPlot(liver_subset, features = B_markers)

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(liver_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Liver cancer Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)


################ 5-Gallbladder cancer subset

# Assign cell types to clusters
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(gallbladder_subset)

gallbladder_subset <- RenameIdents(gallbladder_subset, selected_celltypes)
DimPlot(gallbladder_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(gallbladder_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(gallbladder_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(gallbladder_subset, dims = 1, cells = 500, balanced = TRUE)

# 4. Create Dot Plot
dot_plot <- DotPlot(object = gallbladder_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
ridge_plots <- RidgePlot(gallbladder_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster
VlnPlot(gallbladder_subset, features = B_markers)

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(gallbladder_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Gallbaladder cancer Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)




############### 6 Neuroblastoma

# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(nb_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(nb_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(nb_subset, dims = 1, cells = 500, balanced = TRUE)

# Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(nb_subset)
pbmc_subset <- RenameIdents(pbmc_subset, selected_celltypes)
DimPlot(pbmc_subset, reduction = "umap", label = TRUE, pt.size = 0.5)


# 4. Create Dot Plot
dot_plot <- DotPlot(object = nb_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
ridge_plots <- RidgePlot(nb_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster
VlnPlot(nb_subset, features = B_markers)
# FCRL4 IL2 TCF21

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(nb_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Neuroblastoma Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





########## 7-PBMC (Peripheral Blood Mononuclear Cells
# Define markers of interest
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
DimPlot(pbmc_subset, reduction = "pca") + NoLegend()

# 2. Plot UMAP
DimPlot(pbmc_subset, reduction = "umap")

# 3. Plot Heatmap
DimHeatmap(pbmc_subset, dims = 1, cells = 500, balanced = TRUE)

# Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(pbmc_subset)

pbmc_subset <- RenameIdents(pbmc_subset, selected_celltypes)
DimPlot(pbmc_subset, reduction = "umap", label = TRUE, pt.size = 0.5)

# 4. Create Dot Plot
dot_plot <- DotPlot(object = pbmc_subset, features = B_markers, scale = TRUE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)

# Display the Dot Plot
print(dot_plot)

# 5. Create Ridge Plot for each marker
library(cowplot)
ridge_plots <- RidgePlot(pbmc_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# Display the adjusted Ridge Plot
print(adjusted_ridge_plots)

# 6. Create Violin Plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc_subset, features = B_markers)

# 7. Create Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(pbmc_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "PMC Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





#################### 8-Esophagus

# Définition des marqueurs d'intérêt
B_markers <- c("AICDA", "ASCL1", "CD8B2", "TCF21", "IL21", "RGS13", "FCRL4", "CXCL13", "FOXP3")

# 1. Plot PCA
pca_plot <- DimPlot(esophagus_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(esophagus_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(esophagus_subset, dims = 1, cells = 500, balanced = TRUE)

# 4. Assigner les types de cellules aux clusters après les plots PCA et UMAP
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(esophagus_subset)

esophagus_subset <- RenameIdents(esophagus_subset, selected_celltypes)
umap_with_labels <- DimPlot(esophagus_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)


# 5. Dot Plot
dot_plot <- DotPlot(object = esophagus_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot pour chaque marqueur
ridge_plots <- RidgePlot(esophagus_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
plot(ridge_plots)

# Ajuster les dimensions des plots en utilisant cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")

# 7. Violin Plot - Visualiser les distributions d'expression de chaque cellule unique dans chaque cluster
violin_plot <- VlnPlot(esophagus_subset, features = B_markers)
plot(violin_plot)

# 8. Feature Plot - Visualiser l'expression des caractéristiques dans l'espace de basse dimension
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(esophagus_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Esophagus Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





###################### 9. HNCs


# 1. Plot PCA
pca_plot <- DimPlot(hnsc_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(hnsc_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(hnsc_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(hnsc_subset)

hnsc_subset <- RenameIdents(hnsc_subset, selected_celltypes)
umap_with_labels <- DimPlot(hnsc_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = hnsc_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(hnsc_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(hnsc_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(hnsc_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "HNSC Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





################### 10. Kidney 
# 1. Plot PCA
pca_plot <- DimPlot(kidney_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(kidney_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(kidney_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(kidney_subset)

kidney_subset <- RenameIdents(kidney_subset, selected_celltypes)
umap_with_labels <- DimPlot(kidney_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = kidney_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(kidney_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(kidney_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(kidney_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Kidney Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





############### 11. Skin cancer
# 1. Plot PCA
pca_plot <- DimPlot(skin_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(skin_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(skin_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(skin_subset)

skin_subset <- RenameIdents(skin_subset, selected_celltypes)
umap_with_labels <- DimPlot(skin_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = skin_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(skin_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(skin_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(skin_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = " Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





###################### 12. Lung cancer subset
# 1. Plot PCA
pca_plot <- DimPlot(lung_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(lung_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(lung_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(lung_subset)

lung_subset <- RenameIdents(lung_subset, selected_celltypes)
umap_with_labels <- DimPlot(lung_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = lung_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(lung_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(lung_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(lung_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Lung Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





################ 13. endometrioid_subset

# 1. Plot PCA
pca_plot <- DimPlot(endometrioid_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(endometrioid_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(endometrioid_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(endometrioid_subset)

endometrioid_subset <- RenameIdents(endometrioid_subset, selected_celltypes)
umap_with_labels <- DimPlot(endometrioid_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = endometrioid_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(endometrioid_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(endometrioid_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(endometrioid_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Endometrioid Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





######################## 14.Pancreas
# 1. Plot PCA
pca_plot <- DimPlot(pancreas_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(pancreas_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(pancreas_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(pancreas_subset)

pancreas_subset <- RenameIdents(pancreas_subset, selected_celltypes)
umap_with_labels <- DimPlot(pancreas_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = pancreas_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(pancreas_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(pancreas_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(pancreas_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Pancreas Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)




######################### 15. Thyroid
# 1. Plot PCA
pca_plot <- DimPlot(thyroid_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(thyroid_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(thyroid_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(thyroid_subset)

thyroid_subset <- RenameIdents(thyroid_subset, selected_celltypes)
umap_with_labels <- DimPlot(thyroid_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = thyroid_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(thyroid_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(thyroid_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(thyroid_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Thyroid Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)






############################ 16. Bladder cancer
# 1. Plot PCA
pca_plot <- DimPlot(bladder_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(bladder_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(bladder_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(bladder_subset)

bladder_subset <- RenameIdents(bladder_subset, selected_celltypes)
umap_with_labels <- DimPlot(bladder_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = bladder_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(bladder_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(bladder_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(bladder_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Bladder Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)




######################## 17. Cervix Cancer

# 1. Plot PCA
pca_plot <- DimPlot(cervix_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(cervix_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(cervix_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(cervix_subset)

cervix_subset <- RenameIdents(cervix_subset, selected_celltypes)
umap_with_labels <- DimPlot(cervix_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = cervix_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(cervix_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(cervix_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(cervix_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Cervix Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)






####################### 18. Ovary

# 1. Plot PCA
pca_plot <- DimPlot(ovary_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(ovary_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(ovary_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(ovary_subset)

ovary_subset <- RenameIdents(ovary_subset, selected_celltypes)
umap_with_labels <- DimPlot(ovary_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = ovary_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(ovary_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(ovary_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(ovary_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Ovary Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





######################## 19. Stomach

# 1. Plot PCA
pca_plot <- DimPlot(stomach_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(stomach_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(stomach_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(stomach_subset)

stomach_subset <- RenameIdents(stomach_subset, selected_celltypes)
umap_with_labels <- DimPlot(stomach_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = stomach_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(stomach_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(stomach_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(stomach_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "Stomach Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)





######################### 20. GIST : Gastrointestinal Stromal Tumor
# 1. Plot PCA
pca_plot <- DimPlot(gist_subset, reduction = "pca") + NoLegend()
print(pca_plot)

# 2. Plot UMAP
umap_plot <- DimPlot(gist_subset, reduction = "umap")
print(umap_plot)

# 3. Heatmap
heatmap_plot <- DimHeatmap(gist_subset, dims = 1, cells = 500, balanced = TRUE)
print(heatmap_plot)

# 4. Assign cell types to clusters after PCA and UMAP plots
selected_celltypes <- c("B.01.TCL1A+naiveB", "B.08.ITGB1+SwBm", "B.09.DUSP4+AtM")
names(selected_celltypes) <- levels(gist_subset)

gist_subset <- RenameIdents(gist_subset, selected_celltypes)
umap_with_labels <- DimPlot(gist_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
print(umap_with_labels)

# 5. Dot Plot
dot_plot <- DotPlot(object = gist_subset, features = B_markers, scale = FALSE, group.by = "celltype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = B_markers, labels = B_markers)
print(dot_plot)

# 6. Ridge Plot for each marker
ridge_plots <- RidgePlot(gist_subset, features = B_markers, ncol = 2) + 
  theme(axis.text.y = element_text(size = 4),
        plot.margin = margin(5, 5, 5, 20))
print(ridge_plots)

# Adjust plot dimensions using cowplot
adjusted_ridge_plots <- plot_grid(ridge_plots, ncol = 1, align = "v", axis = "lr")
print(adjusted_ridge_plots)

# 7. Violin Plot - Visualize single cell expression distributions in each cluster
violin_plot <- VlnPlot(gist_subset, features = B_markers)
print(violin_plot)

# 8. Feature Plot - Visualize feature expression in low-dimensional space
# Define a minimum cutoff value to exclude cells with zero expression
min_cutoff_value <- 0.1  # Adjust this value as needed

# FeaturePlots for each marker with adjusted min.cutoff and scale
feature_plots <- lapply(B_markers, function(marker) {
  FeaturePlot(gist_subset, features = marker, min.cutoff = min_cutoff_value, pt.size = 0.5) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "red"), limits = c(0, 3)) +
    ggtitle(paste("Feature plot of", marker)) +
    theme(plot.title = element_text(size = 8),
          axis.title.y = element_text(size = 5),
          axis.title.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "right")
})

# Combine plots using patchwork
combined_plots <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "GIST Feature plots of the 9 genes of interest")

# Display the combined plots
print(combined_plots)



