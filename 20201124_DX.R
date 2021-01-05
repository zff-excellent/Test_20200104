

# Load libraries
.libPaths()
library(easypackages)
packages <- c('Seurat', 'monocle', 'plyr', 'dplyr', 'tidyverse', 'cowplot', 'ggplot2', 'patchwork')
libraries(packages)


###############################################################################################################################################
## Prepare input files

# source("G:/Program_zff/05_Single_cells_2020/Seurat_reanalysis/01_Single_cells_scripts/Standardized_analysis/single_sample/functions.R")
source("I:/Single_Cell_Workflow/Standardized_analysis/functions_2.R")

dataset_path <- 'H:/Cellranger_Data_2020/2020-10-26_N_IRd5A_IRd5cA' # one sample
output_path <- 'G:/Program_zff/05_Single_cells_2020/Seurat_reanalysis/HUST_tongji_Professor_Zeng_Rui/DX/2020_11_06_analysis'
species_names <- 'Mouse' # 'Human' or 'Mouse'.
species_name <- 'mm' # hs or mm


path <- output_path
folder_names <- c('Big_Rdata', '01_Assess_Quality', '02_Clustering_Cells', '03_Differentially_expressed_genes', '04_Celltype_Identification', '05_Compare_between_samples')

path1 <- file.path(path, folder_names[2])
subfolder_names <- c('Raw_data', 'Clean_data')

path2 <- file.path(path, folder_names[5])
subfolder_names2 <- c('SingleR','clustermole')

##### Load seruat object

setwd(file.path(path, folder_names[1]))
# save(data.filt, file='Cluster_all.Rdata', compress="bzip2", compression_level=9)
load("Cluster_all.Rdata")

# step6: Visualize
sapply(grep("res",sort(colnames(data.filt@meta.data)),value = TRUE),function(x) length(unique(data.filt@meta.data[,x])))

Idents(data.filt) <- "integrated_snn_res.0.2"
Idents(data.filt) <- factor(Idents(data.filt), levels = sort(as.numeric(levels(Idents(data.filt)))))
table(Idents(data.filt))


##### Load marker gene list

DefaultAssay(data.filt) <- "RNA"
setwd(file.path(path, folder_names[4]))
markers_all <- read.table("markers_all_DEGs_among_clusters.csv", sep=",", header=T)
top5 <- markers_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) #20
# write.table(top5,"top5_DEGs_among_clusters.csv", sep=",", row.names=F)


###############################################################################################################################################
## clusterProfiler analysis


## Convert geneID
library(clusterProfiler)
library(org.Mm.eg.db)
organism <- 'org.Mm.eg.db'
org <- 'mmu'

gene.df <- bitr(markers_all$gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的  ENTREZID
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = organism)

head(gene.df)


## Load database
# install.packages("vroom")
library(vroom)

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/download.jsp
## All cell markers/Human cell markers/Mouse cell markers/Single cell markers
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers


## Enrichment
y <- enricher(gene.df$ENTREZID, TERM2GENE=cell_markers, minGSSize=1)
DT::datatable(as.data.frame(y))

id <- gene.df$ENTREZID
names(id) <- gene.df$SYMBOL  # 用R语言的names函数构建类似字典的数据结构

## Add gene symbol
y@result-> res  # 以防数据被搞坏，新建一个
res$genesym <- unlist(lapply(y@result$geneID,  FUN =function(x){paste(unlist(lapply(unlist(str_split(x,"/")),FUN=function(x){names(id[which(id ==x)])})) , collapse = "/")} ))
DT::datatable(res)


## Batch for each cluster

for(i in names(table(markers_all$cluster))){
  
  i <- names(table(markers_all$cluster))[1]
  print(i)
  
  ## Load database
  
  ## DEGs
  marker <- FindMarkers(data.filt, ident.1 = i, logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)
  marker <- marker %>% tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, pct.1, pct.2, avg_logFC, p_val, p_val_adj) %>% # FindMarkers
    # mutate_if(is.numeric, round, digits=3) %>%
    dplyr::arrange(-avg_logFC) %>% 
    subset(p_val_adj < 0.01)
  print(dim(marker))
  write.table(marker, file=paste0("C", i, "_marker_genes.csv"), sep=",", row.names=F)
  
  ## Convert geneID
  genes <- subset(markers_all, cluster==i) %>% top_n(100, wt=avg_logFC)
  gene.df <- bitr(genes$gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的  ENTREZID
                  toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                  OrgDb = organism)
  print(dim(gene.df)) ## 99
  
  ## Enrichment
  y <- enricher(gene.df$ENTREZID, TERM2GENE=cell_markers, minGSSize=1) 
  DT::datatable(as.data.frame(y))
  
  ## Add gene symbol
  y@result -> res  # 以防数据被搞坏，新建一个
  res$genesym <- unlist(lapply(y@result$geneID,  FUN =function(x){paste(unlist(lapply(unlist(str_split(x,"/")),FUN=function(x){names(id[which(id ==x)])})) , collapse = "/")} ))
  DT::datatable(res)
  
  res$geneID
  
  ## Save files
  write.csv(as.data.frame(res[-1]), file=paste0("C", i, "_KEGG_annotation.csv"),row.names =F)
  
}


cd_genes <- unlist(str_split(strings, '/'))
unlist(str_split(res$geneID, '/'))
unique(sort(unlist(str_split(res$geneID, '/'))))



