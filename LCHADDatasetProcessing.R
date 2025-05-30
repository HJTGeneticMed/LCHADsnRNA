####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Dataset Initial Processing, QC, Filtering Cells, Doublet Identification, and scVI Integration #### Pipeline used https://github.com/XL-Genomics/2022_macrophage_cardiac_allograft/tree/main STEP 01-07
####  Go to: https://raw.githubusercontent.com/XL-Genomics/2022_macrophage_cardiac_allograft/main/src/bioinformatics.R, https://github.com/XL-Genomics/2022_macrophage_cardiac_allograft/blob/main/src/scRNAseq.R, and run the scripts in R. 
####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Global Annotation (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHAD <- RunUMAP(LCHAD, reduction = "scVI2", dims = 1:50)
LCHAD <- FindNeighbors(LCHAD reduction = "scVI2", dims = 1:50)
LCHAD <- FindClusters(LCHAD, resolution = c(0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHAD, reduction = 'umap', group.by = 'Doublet_SC')

DimPlot2(LCHAD, reduction = 'umap', group.by = 'RNA_snn_res.0.1', label = TRUE) +
  DimPlot2(LCHAD, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHAD, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE)

p1 <- DotPlot2(LCHAD, features = grep('^Markerset_', colnames(LCHAD@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHAD, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHAD$Cell_type <- NA
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(1, 12, 14, 17)] <- "CM"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(0, 5, 21)]       <- "FB"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(2, 6, 7, 19)]    <- "Endo"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(3)]              <- "Peri"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(4, 23, 24)]      <- "Myeloid"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(8)]              <- "Lymph"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(9)]              <- "Neuro"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(10)]             <- "SMC"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(18)]             <- "Meso"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(26)]             <- "Adipocyte"
LCHAD$Cell_type[LCHAD$RNA_snn_res.0.5 %in% c(11, 13, 15, 16, 20, 22, 25)] <- "Doublet"

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHAD, 'LCHADPreSubClustering.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD CM Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADCM <- LCHAD[, LCHAD$Cell_type %in% c('CM')]
LCHADCM <- RunUMAP(LCHADCM, reduction = "scVI2", dims = 1:50)
LCHADCM <- FindNeighbors(LCHADCM reduction = "scVI2", dims = 1:50)
LCHADCM <- FindClusters(LCHADCM, resolution = c(0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADCM, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADCM, reduction = umap', group.by = 'RNA_snn_res.0.1', label = T) +
        DimPlot2(LCHADCM, reduction = ‘umap', group.by = 'RNA_snn_res.0.2', label = T) +
  DimPlot2(LCHADCM, reduction = ' umap', group.by = 'RNA_snn_res.0.5', label = T)
p1 <- DotPlot2(LCHADCM, features = grep('^Markerset_', colnames(LCHADCM@meta.data), value = T))

Markersforfirstscreen <- FindAllMarkers(LCHADCM, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHADCM$Cell_state <- NA
LCHADCM$Cell_state <- as.vector(LCHADCM$RNA_snn_res.0.5)
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(0)] <- 'CM1'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(1)] <- 'CM2'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(2)] <- 'CM3'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(3)] <- 'CM4'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(5)] <- 'CM5'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(6)] <- 'CM6'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(7)] <- 'CM7'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(8)] <- 'CM8'
LCHADCM$Cell_state[LCHADCM$RNA_snn_res.0.5 %in% c(4, 9, 10)] <- 'Doublet'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADCM, 'LCHADCM.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD FB Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADFB <- LCHAD[, LCHAD$Cell_type %in% c('FB')] 
LCHADFB <- RunUMAP(LCHADFB, reduction = "scVI2", dims = 1:50)
LCHADFB <- FindNeighbors(LCHADFB, reduction = "scVI2", dims = 1:50)
LCHADFB <- FindClusters(LCHADFB, resolution = c(0.2, 0.3, 0.35, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADFB, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADFB, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHADFB, reduction = 'umap', group.by = 'RNA_snn_res.0.3', label = TRUE) +
  DimPlot2(LCHADFB, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE)
p1 <- DotPlot2(LCHADFB, features = grep('^Markerset_', colnames(LCHADFB@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHADFB, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
Idents(LCHADFB) <- "RNA_snn_res.0.35"
LCHADFB$Cell_state <- as.vector(LCHADFB$RNA_snn_res.0.35)
LCHADFB$Cell_state[LCHADFB$RNA_snn_res.0.35 %in% c(0)] <- 'FB2'
LCHADFB$Cell_state[LCHADFB$RNA_snn_res.0.35 %in% c(1)] <- 'FB4'
LCHADFB$Cell_state[LCHADFB$RNA_snn_res.0.35 %in% c(2)] <- 'FB1'
LCHADFB$Cell_state[LCHADFB$RNA_snn_res.0.35 %in% c(3)] <- 'FB3'
LCHADFB$Cell_state[LCHADFB$RNA_snn_res.0.35 %in% c(4, 5, 6)] <- 'Doublet'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADFB, 'LCHADFB.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD ECs Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADEC <- LCHAD[, LCHAD$Cell_type %in% c('Endo')] 
LCHADEC <- RunUMAP(LCHADEC, reduction = "scVI2", dims = 1:50)
LCHADEC <- FindNeighbors(LCHADEC, reduction = "scVI2", dims = 1:50)
LCHADEC <- FindClusters(LCHADEC, resolution = c(0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADEC, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADEC, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHADEC, reduction = 'umap', group.by = 'RNA_snn_res.0.3', label = TRUE) +
  DimPlot2(LCHADEC, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE)
p1 <- DotPlot2(LCHADEC, features = grep('^Markerset_', colnames(LCHADEC@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHADEC, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHADEC$Cell_state <- NA
LCHADEC$Cell_state[LCHADEC$RNA_snn_res.0.2 %in% 0] <- 'EC1'
LCHADEC$Cell_state[LCHADEC$RNA_snn_res.0.2 %in% 1] <- 'EC2'
LCHADEC$Cell_state[LCHADEC$RNA_snn_res.0.2 %in% 2] <- 'EC3'
LCHADEC$Cell_state[LCHADEC$RNA_snn_res.0.2 %in% 3] <- 'Endo'
LCHADEC$Cell_state[LCHADEC$RNA_snn_res.0.2 %in% 4] <- 'LECs'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADEC, 'LCHADEC.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD ML/Lymph Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADMLLymph <- LCHAD[, LCHAD$Cell_type %in% c('Myeloid', 'Lymph')] 
LCHADMLLymph <- RunUMAP(LCHADMLLymph, reduction = "scVI2", dims = 1:50)
LCHADMLLymph <- FindNeighbors(LCHADMLLymph, reduction = "scVI2", dims = 1:50)
LCHADMLLymph <- FindClusters(LCHADMLLymph, resolution = c(0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADMLLymph, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADMLLymph, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHADMLLymph, reduction = 'umap', group.by = 'RNA_snn_res.0.3', label = TRUE) +
  DimPlot2(LCHADMLLymph, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE)
p1 <- DotPlot2(LCHADMLLymph, features = grep('^Markerset_', colnames(LCHADMLLymph@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHADMLLymph, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHADMLLymph$Cell_state <- NA
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(1)] <- "MP1"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(0)] <- "MP2"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(6)] <- "Mo-HLAhi"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(3)] <- "LAM-like"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(9)] <- "DC"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(7)] <- "Mast"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(2)] <- "CD4+ TC"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(5)] <- "NK"
LCHADMLLymph$Cell_state[LCHADMLLymph$RNA_snn_res.0.5 %in% c(4, 8)] <- "Doublet"

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADMLLymph, 'LCHADMLLymph.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD SMC Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADSMC <- LCHAD[, LCHAD$Cell_type %in% c('SMC')] 
LCHADSMC <- RunUMAP(LCHADSMC, reduction = "scVI2", dims = 1:50)
LCHADSMC <- FindNeighbors(LCHADSMC, reduction = "scVI2", dims = 1:50)
LCHADSMC <- FindClusters(LCHADSMC, resolution = c(0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADSMC, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADSMC, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHADSMC, reduction = 'umap', group.by = 'RNA_snn_res.0.3', label = TRUE) +
  DimPlot2(LCHADSMC, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE)
p1 <- DotPlot2(LCHADSMC, features = grep('^Markerset_', colnames(LCHADSMC@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHADSMC, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHADSMC$Cell_state <- NA
LCHADSMC$Cell_state[LCHADSMC$RNA_snn_res.0.2 %in% c(0)] <- "SMC1"
LCHADSMC$Cell_state[LCHADSMC$RNA_snn_res.0.2 %in% c(1)] <- "SMC2"

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADSMC, 'LCHADSMC.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Peri Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADPC <- LCHAD[, LCHAD$Cell_type %in% c('Peri')] 
LCHADPC <- RunUMAP(LCHADPC, reduction = "scVI2", dims = 1:50)
LCHADPC <- FindNeighbors(LCHADPC, reduction = "scVI2", dims = 1:50)
LCHADPC <- FindClusters(LCHADPC, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

## Cluster Doublet Identification
DimPlot2(LCHADPC, reduction = 'umap', group.by = 'Doublet_SC')
DimPlot2(LCHADPC, reduction = 'umap', group.by = 'RNA_snn_res.0.2', label = TRUE) +
  DimPlot2(LCHADPC, reduction = 'umap', group.by = 'RNA_snn_res.0.3', label = TRUE) +
  DimPlot2(LCHADPC, reduction = 'umap', group.by = 'RNA_snn_res.0.1', label = TRUE)
p1 <- DotPlot2(LCHADPC, features = grep('^Markerset_', colnames(LCHADSMC@meta.data), value = TRUE))

Markersforfirstscreen <- FindAllMarkers(LCHADPC, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
Markersforfirstscreen <- Markersforfirstscreen[Markersforfirstscreen$p_val_adj < 0.05, ]

## Cluster Re-Annotation
LCHADPC$Cell_state <- NA
LCHADPC$Cell_state[LCHADPC$RNA_snn_res.0.1 %in% c(0)] <- "PC1"
LCHADPC$Cell_state[LCHADPC$RNA_snn_res.0.1 %in% c(1)] <- "PC2"

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADPC, 'LCHADPC.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Neuro Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADNeuro <- LCHAD[, LCHAD$Cell_type %in% c('Neuro')] 
LCHADNeuro <- RunUMAP(LCHADNeuro, reduction = "scVI2", dims = 1:50)
LCHADNeuro <- FindNeighbors(LCHADNeuro, reduction = "scVI2", dims = 1:50)
LCHADNeuro <- FindClusters(LCHADNeuro, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

## Cluster Re-Annotation; small nuclei number
LCHADNeuro$Cell_state <- 'Neuro'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADNeuro, 'LCHADNeuro.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Adipo Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADAdipo <- LCHAD[, LCHAD$Cell_type %in% c('Adipo')] 
LCHADAdipo <- RunUMAP(LCHADAdipo, reduction = "scVI2", dims = 1:50)
LCHADAdipo <- FindNeighbors(LCHADAdipo, reduction = "scVI2", dims = 1:50)
LCHADAdipo <- FindClusters(LCHADAdipo, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

## Cluster Re-Annotation; small nuclei number
LCHADAdipo$Cell_state <- 'Adipo'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADAdipo, 'LCHADAdipo.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD Mesothelial Cells Subclustering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

## Re-embed by scVI
LCHADMeso <- LCHAD[, LCHAD$Cell_type %in% c('Meso')] 
LCHADMeso <- RunUMAP(LCHADMeso, reduction = "scVI2", dims = 1:50)
LCHADMeso <- FindNeighbors(LCHADMeso, reduction = "scVI2", dims = 1:50)
LCHADMeso <- FindClusters(LCHADMeso, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

## Cluster Re-Annotation; small nuclei number
LCHADMeso$Cell_state <- 'Meso'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHADMeso, 'LCHADMeso.rds')
####--------------------------------------------------------------------------------------------------------------------

####--------------------------------------------------------------------------------------------------------------------
####  Final LCHAD Global Processing following SubClutering (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------

LCHAD <- readRDS("~/LCHAD/srt.rds")

LCHAD$Cell_state <- NA
LCHAD$Cell_type <- as.vector(LCHAD$Cell_type)

LCHAD@reductions$sub_umap@cell.embeddings[,1] <- NA
LCHAD@reductions$sub_umap@cell.embeddings[,2] <- NA
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADCM),] <- LCHADCM@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADFB),] <- LCHADFB@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADEC),] <- LCHADEC@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADMLLymph),] <- LCHADMLLymph@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADPC),] <- LCHADPC@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADSMC),] <- LCHADSMC@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADMeso),] <- LCHADMeso@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADNeuro),] <- LCHADNeuro@reductions$umap@cell.embeddings
LCHAD@reductions$sub_umap@cell.embeddings[Cells(LCHADAdipo),] <- LCHADAdipo@reductions$umap@cell.embeddings

LCHAD$Cell_state[Cells(LCHADCM)] <- as.vector(LCHADCM$Cell_state)
LCHAD$Cell_state[Cells(LCHADFB)] <- as.vector(LCHADFB$Cell_state)
LCHAD$Cell_state[Cells(LCHADEC)] <- as.vector(LCHADEC$Cell_state)
LCHAD$Cell_state[Cells(LCHADMLLymph)] <- as.vector(LCHADMLLymph$Cell_state)
LCHAD$Cell_state[Cells(LCHADPC)] <- as.vector(LCHADPC$Cell_state)
LCHAD$Cell_state[Cells(LCHADSMC)] <- as.vector(LCHADSMC$Cell_state)
LCHAD$Cell_state[Cells(LCHADMeso)] <- as.vector(LCHADMeso$Cell_state)
LCHAD$Cell_state[Cells(LCHADNeuro)] <- as.vector(LCHADNeuro$Cell_state)
LCHAD$Cell_state[Cells(LCHADAdipo)] <- as.vector(LCHADAdipo$Cell_state)

LCHAD$Cell_type[LCHAD$Cell_state %in% c('CD4 TC', 'NK')] <- 'Lymph'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('CM1', 'CM2', 'CM3', 'CM4', 'CM5', 'CM6', 'CM7', 'CM8')] <- 'CM'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('FB1', 'FB2', 'FB3', 'FB4')] <- 'FB'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('EC1', 'EC2', 'EC3', 'Endo', 'LECs')] <- 'EC'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('MP1', 'MP2', 'Mo-HLAhi', 'LAM-like', 'DC', 'Mast')] <- 'ML'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('PC1', 'PC2')] <- 'PC'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('SMC1', 'SMC2')] <- 'SMC'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('Meso')] <- 'Meso'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('Neuro')] <- 'Neuro'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('Adipo')] <- 'Adipo'
LCHAD$Cell_type[LCHAD$Cell_state %in% c('Doublet')] <- 'Doublet'

####--------------------------------------------------------------------------------------------------------------------
saveRDS(LCHAD, 'LCHADProcessed.rds')
####--------------------------------------------------------------------------------------------------------------------