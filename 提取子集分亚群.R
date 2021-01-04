
library(easypackages)
packages <- c('ggplot2', 'cowplot', 'plyr', 'dplyr', 'tidyverse', 'patchwork', 'pheatmap','Seurat', 'monocle', 'scRNAseq', 
              'SingleR',  'scCATCH', 'scMCA', 'scHCL', 'clustermole')
# install.packages(packages)
libraries(packages)


setwd('G:/Program_zff/05_Single_cells_2020/Seurat_reanalysis/Xiehe_hospital_Professor_Puhongxu/2020_12_07_IFP1-SCAT1')
load('./Rdata/Cluster_all.Rdata')

Idents(data.filt) <- "integrated_snn_res.0.8"
Idents(data.filt) <- factor(Idents(data.filt), levels = sort(as.numeric(levels(Idents(data.filt)))))
table(Idents(data.filt))

###############################################################################################################################################
## Rename celltype and visualize genes distribution

label_names <- read.csv('Rename_Label.csv', header=F)
labers = label_names[match(as.numeric(as.character(data.filt@active.ident)),label_names[,1]),2]
data.filt$labers <- labers
# DimPlot(object = data.filt, reduction = "tsne", group.by = "labers", label = TRUE, pt.size=1)
p1 <- DimPlot(object = data.filt, reduction = "tsne", group.by = "labers", label = TRUE, pt.size=1) 
p2 <- DimPlot(object = data.filt, reduction = "umap", group.by = "labers", label = TRUE, pt.size=1)
p <- plot_grid(p1, p2)
p

table(Idents(data.filt))
table(data.filt$labers)


###########################################################
## Way

DefaultAssay(data.filt) <- 'RNA'
sub_cells <- subset(data.filt, labers == 'ADSC')

# DefaultAssay(sub_cells) <- 'RNA'
pd <- sub_cells@meta.data
new_pd <- dplyr::select(pd,orig.ident,nCount_RNA,nFeature_RNA,percent.mito)
sub_cells@meta.data <- new_pd

# step1: Normalization
seurat_list <- SplitObject(sub_cells, split.by = "orig.ident")

for (i in 1:length(seurat_list)) {
  # sparse_data <- as(seurat_list[[i]][["RNA"]]@data, "dgCMatrix") # seurat to dgcmatrix, then filter
  # seurat_list[[i]] <- CreateSeuratObject(counts = sparse_data, min.cells = 3, min.features = 200)
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000) # from 2000 to 3000
} 

# step2, Integrate samples using shared highly variable genes

# Identify anchors
# # This appear to be related to the number of cells in the smallest dataset. 
# scseqs <- seurat_list
# k.filter <- min(200, min(sapply(scseqs, ncol)))
# anchors <- Seurat::FindIntegrationAnchors(scseqs, k.filter = k.filter)

anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)

# Integrate samples
combined <- IntegrateData(anchorset = anchors)

length(VariableFeatures(combined))

# step3: normalization together
sub_cells <- combined
sub_cells <- ScaleData(sub_cells, verbose = FALSE)# features = rownames(combined), vars.to.regress = c("percent.mito", "percent.ribo", "orig.ident")
sub_cells <- RunPCA(sub_cells, verbose = FALSE) # features = rownames(combined)

seurat_control <- sub_cells
pct <- seurat_control[["pca"]]@stdev / sum(seurat_control[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- max(co1, co2)
use.pcs = 1:pcs

# step4: Determine the K-nearest neighbor graph
sub_cells <- FindNeighbors(sub_cells, reduction="pca", dims = use.pcs)
# Determine the clusters for various resolutions
sub_cells <- FindClusters(sub_cells, resolution = seq(0.1,5,0.1))

# step5: t-SNE and Clustering
sub_cells <- RunTSNE(sub_cells, reduction = "pca", dims = use.pcs) # check_duplicates = FALSE, Remove duplicates before running TSNE
sub_cells <- RunUMAP(sub_cells, reduction = "pca", dims = use.pcs)


save(sub_cells, file = 'Sub_Cluster_all.Rdata')


sapply(grep("res",colnames(sub_cells@meta.data),value = TRUE),function(x) length(unique(sub_cells@meta.data[,x])))
Idents(sub_cells) <- 'integrated_snn_res.0.2'
table(Idents(sub_cells))

p1 <- DimPlot(object = sub_cells, reduction = "tsne", label = TRUE, pt.size=1) 
p2 <- DimPlot(object = sub_cells, reduction = "umap", label = TRUE, pt.size=1)
p <- plot_grid(p1, p2)
p

save_plot("Compare_TSNE_and_UMAP.png", p, base_height = 8, base_aspect_ratio = 2.5, base_width = NULL, dpi=600) # 2.5/3.5
save_plot("Compare_TSNE_and_UMAP.pdf", p, base_height = 8, base_aspect_ratio = 2.5, base_width = NULL)

save_plot("TSNE_clustering.png", p1, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL, dpi=600)
save_plot("TSNE_clustering.pdf", p1, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL)
save_plot("UMAP_clustering.png", p2, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL, dpi=600)
save_plot("UMAP_clustering.pdf", p2, base_height = 8, base_aspect_ratio = 1.3, base_width = NULL)


# Visualize marker genes
markers_all <- FindAllMarkers(object = sub_cells, min.pct = 0.1, logfc.threshold = 0.25, only.pos=T) # only.pos=T, min.pct=0.25
markers_all <- dplyr::select(markers_all, gene, cluster, pct.1, pct.2, avg_logFC, p_val, p_val_adj)
write.table(markers_all,"markers_all_DEGs_among_clusters.csv", sep=",", row.names=F)


###########################################################
###########################################################
