####--------------------------------------------------------------------------------------------------------------------
####  Downstream Analyses
####  In order to run functions, go to: https://raw.githubusercontent.com/XL-Genomics/2022_macrophage_cardiac_allograft/main/src/bioinformatics.R, https://github.com/XL-Genomics/2022_macrophage_cardiac_allograft/blob/main/src/scRNAseq.R, and run the scripts in R. 
####--------------------------------------------------------------------------------------------------------------------
####  LCHAD CM Downstream Analyses ####
####--------------------------------------------------------------------------------------------------------------------

LCHADCM <- readRDS("~/LCHAD/LCHADCM.rds")

## Cluster Gene Markers
LCHADCM <- LCHADCM[, !(LCHADCM$Cell_state %in% c("Doublet"))]
Idents(LCHADCM) <- LCHADCM$Cell_state
FinalCMMarkers <- FindAllMarkers(LCHADCM, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalCMMarkers <- FinalCMMarkers[FinalCMMarkers$p_val_adj < 0.05, ]
write.csv(FinalCMMarkers, 'FinalLCHADCMMarkers.csv')

## CM Cell State BarPlot
p1 <- CountCellBarPlot(LCHADCM, group.var = 'group2', stack.var = 'Cell_state', percentage = TRUE, stack.color = mycol_10)$plot

p1 <- p1 +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    axis.line.y = element_line(size = 0.8, color = "black"),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 12, color = "black", face = "plain"),
    axis.text.y = element_text(size = 12, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black", face = "plain", margin = margin(t = -4)),
    legend.text = element_text(size = 12, color = "black", face = "plain"),
    legend.title = element_text(size = 12, color = "black", face = "plain"),
    legend.key.size = unit(1.3, "lines"),
    legend.spacing = unit(2, "cm")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_x_discrete(labels = c("Control", "Ped DCM", "Ped LCHAD")) +
  scale_fill_manual(name = "Cell State", values = mycol_10) +
  labs(
    y = "Percentage",
    x = ""
  ) +
  guides(fill = guide_legend(override.aes = list(colour = "black")))

print(p1)

PlotPDF('LCADHCMCellStateBarPlot', 4, 5)
p1
dev.off()

## Pathway Enrichment Analysis by Group

Idents(LCHADCM) <- LCHADCM$group2
FinalCMbyGroup <- FindAllMarkers(LCHADCM, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalCMbyGroup <- FinalCMbyGroup[FinalCMbyGroup$p_val_adj < 0.05, ]
FinalCMbyGroup <- FinalCMbyGroup[FinalCMbyGroup$avg_log2FC > 0.5, ]
write.csv(FinalCMbyGroup, 'FinalCMbyGroupMarkers.csv')
genelist <- split(FinalCMbyGroup$gene, FinalCMbyGroup$cluster)
all_genes <- unique(unlist(genelist))
protein_coding_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = mart)
protein_coding_genes <- protein_coding_info$hgnc_symbol[
  protein_coding_info$gene_biotype == "protein_coding"
]
filtered_genelist <- lapply(genelist, function(glist) {
  glist[glist %in% protein_coding_genes]
})
filtered_genelist
for(name in names(filtered_genelist)) {
  df <- data.frame(gene = filtered_genelist[[name]])
  filename <- paste0(name, "_filtered_genes.csv")
  write.csv(df, file = filename, row.names = FALSE)}

## Module Score Creation

GeneModuleScores <- read_excel("~Kuppe.Nature.GeneModuleScores.xlsx")

LCHADCM <- AddModuleScore2(
  LCHADCM,
  features = list(GeneModuleScores$`Stressed Myocardium Module Score Genes`),
  name = "Stressed",
  return_z = TRUE
)

## Normal and Stress Myocardial Scores
LCHADCM <- AddModuleScore2(
  LCHADCM,
  features = list(GeneModuleScores$`Normal Myocardium Module Score Genes`),
  name = "Normal",
  return_z = TRUE)

p1 <- VlnPlot2(
  LCHADCM,
  features = "Normal",
  group.by = "group2",
  pt.size = 0, col = mycol_10) +
  geom_boxplot(
    width = 0.3,
    outlier.shape = NA,
    color = "black",
    position = position_dodge(width = 5)) +
  labs(x = "Cell State", y = "Normal Z Score") +
  theme(
    axis.title.x = element_blank(),           # Remove x-axis title
    axis.text.x = element_blank(),            # Remove x-axis text
    axis.ticks.x = element_blank(),           # Remove x-axis ticks
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    plot.title = element_blank(),
    legend.position = "none")
# Plot for Stressed
p2 <- VlnPlot2(
  LCHADCM,
  features = "Stressed",
  group.by = "group2",
  pt.size = 0, col = mycol_10) +
  geom_boxplot(
    width = 0.3,
    outlier.shape = NA,
    color = "black",
    position = position_dodge(width = 5)) +
  labs(x = "Cell State", y = "Stressed Z Score") +
  scale_x_discrete(labels = c("Peds_Ctrl" = "Control", "Peds_DCM" = "Ped DCM", "Peds_LCHAD" = "Ped LCHAD")) +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_blank(),
    legend.position = "none")
# Combine
a <- p1 / p2
a
PlotPDF('StressNormalScoreVlnPlot', 3, 6)
a
dev.off()

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD FB Downstream Analyses ####
####--------------------------------------------------------------------------------------------------------------------

LCHADFB <- readRDS("~/LCHAD/LCHADFB.rds")

## Cluster Gene Markers
LCHADFB <- LCHADFB[, !(LCHADFB$Cell_state %in% c("Doublet"))]
FinalFBMarkers <- FindAllMarkers(LCHADFB, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalFBMarkers <- FinalFBMarkers[FinalFBMarkers$p_val_adj < 0.05, ]
write.csv(FinalFBMarkers, "FinalLCHADFBMarkers.csv")

## FB Cell State BarPlot
p1 <- CountCellBarPlot(LCHADFB, group.var = 'group2', stack.var = 'Cell_state', percentage = TRUE, stack.color = mycol_10)$plot

p1 <- p1 +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    axis.line.y = element_line(size = 0.8, color = "black"),
    axis.line.x = element_blank(),
    axis.title = element_text(size = 12, color = "black", face = "plain"),
    axis.text.y = element_text(size = 12, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black", face = "plain", margin = margin(t = -4)),
    legend.text = element_text(size = 12, color = "black", face = "plain"),
    legend.title = element_text(size = 12, color = "black", face = "plain"),
    legend.key.size = unit(1.3, "lines"),
    legend.spacing = unit(2, "cm")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_x_discrete(labels = c("Control", "Ped DCM", "Ped LCHAD")) +
  scale_fill_manual(name = "Cell State", values = mycol_10) +
  labs(
    y = "Percentage",
    x = ""
  ) +
  guides(fill = guide_legend(override.aes = list(colour = "black")))

print(p1)

PlotPDF('LCHADFBCellStateBarPlot', 4, 5)
p1
dev.off()

## Pathway Enrichment Analysis by Group
Idents(LCHADFB) <- LCHADFB$group2
FinalFBbyGroup <- FindAllMarkers(LCHADFB, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalFBbyGroup <- FinalFBbyGroup[FinalFBbyGroup$p_val_adj < 0.05, ]
FinalFBbyGroup <- FinalFBbyGroup[FinalFBbyGroup$avg_log2FC > 0.5, ]
write.csv(FinalFBbyGroup, 'FinalFBbyGroupMarkers.csv')
genelist <- split(FinalFBbyGroup$gene, FinalFBbyGroup$cluster)
all_genes <- unique(unlist(genelist))
protein_coding_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = mart)
protein_coding_genes <- protein_coding_info$hgnc_symbol[
  protein_coding_info$gene_biotype == "protein_coding"
]
filtered_genelist <- lapply(genelist, function(glist) {
  glist[glist %in% protein_coding_genes]
})
filtered_genelist
for(name in names(filtered_genelist)) {
  df <- data.frame(gene = filtered_genelist[[name]])
  filename <- paste0(name, "_filtered_genes.csv")
  write.csv(df, file = filename, row.names = FALSE)
}

## Module Score Creation

LCHADFB <- AddModuleScore2(
  LCHADFB,
  features = list(c("COL22A1", "POSTN", "THBS4", "TSHZ2", "FAM155A", "COL1A2", "MIR100HG", "NOX4", "FAP", "COL1A1")),
  name = "ActFBScore",
  return_z = TRUE
)

## ActFB Score
a <- VlnPlot2(LCHADFB, features = c('ActFBScore'), same.y.lims = 6,
              group.by = 'group2', pt.size = 0, ncol = 1, cols = mycol_10) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black",
               position = position_dodge(width = 5)) +
  labs(x = "Cell State", y = "ActFB Score") +  # Add axis titles
  theme(
    axis.title.x = element_text(size = 14),   # Make x-axis title bigger
    axis.title.y = element_text(size = 14),   # Make y-axis title bigger  # Remove x-axis text
    axis.text.y = element_text(size = 14),axis.text.x = element_text(size = 16),    # Make y-axis text bigger
    plot.title = element_blank(),  # Remove plot title
    legend.position = "none"  # Remove the legend entirely
  )
a

PlotPDF('ActFBScoreVlnPlot', 3, 6)
a
dev.off()

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD EC Downstream Analyses ####
####--------------------------------------------------------------------------------------------------------------------

LCHADEC <- readRDS("~/LCHAD/LCHADEC.rds")

## Cluster Gene Markers
LCHADEC <- LCHADEC[, !(LCHADEC$Cell_state %in% c("Doublet"))]
FinalECMarkers <- FindAllMarkers(LCHADEC, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalECMarkers <- FinalECMarkers[FinalECMarkers$p_val_adj < 0.05, ]
write.csv(FinalECMarkers, "FinalECMarkersMarkers.csv")

## Pathway Enrichment Analysis by Cell State
Idents(LCHADEC) <- LCHADEC$group2
FinalECbyGroup <- FindAllMarkers(LCHADEC, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalECbyGroup <- FinalECbyGroup[FinalECbyGroup$p_val_adj < 0.05, ]
FinalECbyGroup <- FinalECbyGroup[FinalECbyGroup$avg_log2FC > 0.5, ]
write.csv(FinalECbyGroup, 'FinalECbyGroupMarkers.csv')
genelist <- split(FinalECbyGroup$gene, FinalECbyGroup$cluster)
all_genes <- unique(unlist(genelist))
protein_coding_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = mart)
protein_coding_genes <- protein_coding_info$hgnc_symbol[
  protein_coding_info$gene_biotype == "protein_coding"
]
filtered_genelist <- lapply(genelist, function(glist) {
  glist[glist %in% protein_coding_genes]
})
filtered_genelist
for(name in names(filtered_genelist)) {
  df <- data.frame(gene = filtered_genelist[[name]])
  filename <- paste0(name, "_filtered_genes.csv")
  write.csv(df, file = filename, row.names = FALSE)
}

## GO Biological Processes Bar Plot
## Input Genes from Mo-HLAhi into https://maayanlab.cloud/Enrichr/ and select for pathways of interest
Top5select <- read_excel("Top5LCHADECselect.xlsx")

top5_pathways <- Top5select %>%
  mutate(`Adjusted P-value` = as.numeric(`Adjusted P-value`),
         neg_log10_adj_pvalue = -log10(`Adjusted P-value`)) %>%
  arrange(desc(neg_log10_adj_pvalue)) %>%
  head(5)
a <- ggplot(top5_pathways, aes(x = neg_log10_adj_pvalue, y = reorder(Term, neg_log10_adj_pvalue))) +
  geom_bar(stat = "identity", fill = "#5AB05AFF") +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_text(aes(x = 0, label = Term), hjust = 0, nudge_x = 0.1, size = 7, color = "black") +
  labs(x = "-Log10(Adjusted P-value)", y = "Enriched GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black")
  )
PlotPDF('ECLCHADGOCC', 5.7, 4.3)
a
dev.off()

####--------------------------------------------------------------------------------------------------------------------
####  LCHAD ML and Lymph Downstream Analyses ####
####--------------------------------------------------------------------------------------------------------------------

LCHADMLLymph <- readRDS("~/LCHAD/LCHADMLLymph.rds")

## Cluster Gene Markers
LCHADMLLymph <- LCHADMLLymph[, !(LCHADMLLymph$Cell_state %in% c("Doublet"))]
FinalMLLymphMarkers <- FindAllMarkers(LCHADMLLymph, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
FinalMLLymphMarkers <- FinalMLLymphMarkers[FinalMLLymphMarkers$p_val_adj < 0.05, ]
write.csv(FinalMLLymphMarkers, "FinalLCHADMLLymphMarkers.csv")

## Pathway Enrichment Analysis by Cell State
Idents(LCHADMLLymph) <- LCHADMLLymph$Cell_state
LCHADMLLymph <- FindAllMarkers(LCHADMLLymph, min.pct = 0.05, logfc.threshold = 0.15, only.pos = TRUE)
LCHADMLLymph <- LCHADMLLymph[LCHADMLLymph$p_val_adj < 0.05, ]
LCHADMLLymph <- LCHADMLLymph[LCHADMLLymph$avg_log2FC > 0.5, ]
genelist <- split(LCHADMLLymph$gene, LCHADMLLymph$cluster)
all_genes <- unique(unlist(genelist))
protein_coding_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = mart)
protein_coding_genes <- protein_coding_info$hgnc_symbol[
  protein_coding_info$gene_biotype == "protein_coding"
]
filtered_genelist <- lapply(genelist, function(glist) {
  glist[glist %in% protein_coding_genes]
})
filtered_genelist
for(name in names(filtered_genelist)) {
  df <- data.frame(gene = filtered_genelist[[name]])
  filename <- paste0(name, "_filtered_genes.csv")
  write.csv(df, file = filename, row.names = FALSE)
}

## UMAP Dimplot
p1 <- DimPlot2(LCHADMLLymph,
               reduction = 'umap.cca',
               group.by = 'Cell_state',
               cols = mycol_10,
               label = FALSE,
               raster = FALSE,
               pt.size = 0.2,
               order = TRUE,
               label.size = 10) +
  labs(title = 'Cell Types', x = 'UMAP1', y = 'UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 6, shape = 21, fill = mycol_10[1:8], color = "black", stroke = 1), ncol = 1)) +
  theme(legend.text = element_text(size = 15))

PlotPDF('LCHADMLLymph.umap', 6, 8)
p1
dev.off()

## GO Biological Processes Bar Plot
## Input Genes from Mo-HLAhi into https://maayanlab.cloud/Enrichr/ and select for pathways of interest
Top5select <- read_excel("Top5MoHLAhiselect.xlsx")
top5_pathways <- Top5select %>%
  mutate(`Adjusted P-value` = as.numeric(`Adjusted P-value`),
         neg_log10_adj_pvalue = -log10(`Adjusted P-value`)) %>%
  arrange(desc(neg_log10_adj_pvalue)) %>%
  head(5)
a <- ggplot(top5_pathways, aes(x = neg_log10_adj_pvalue, y = reorder(Term, neg_log10_adj_pvalue))) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_text(aes(x = 0, label = Term), hjust = 0, nudge_x = 0.1, size = 7, color = "black") +
  labs(x = "-Log10(Adjusted P-value)", y = "Enriched GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black")
  )
PlotPDF('MoHLAhiGOBP', 5.7, 4.3)
a
dev.off()

####--------------------------------------------------------------------------------------------------------------------
####  Complete LCHAD Dataset Figures ####
####--------------------------------------------------------------------------------------------------------------------

LCHAD <- LCHAD[, !(LCHAD$Cell_type %in% c("Doublet"))]

p1 <- DimPlot2(LCHAD,
               reduction = 'umap.cca',
               group.by = 'Cell_type',
               cols = mycol_10,
               label = FALSE,
               raster = FALSE,
               pt.size = 0.2,
               order = TRUE,
               label.size = 10) +
  labs(title = 'Cell Types', x = 'UMAP1', y = 'UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 6, shape = 21, fill = mycol_10[1:10], color = "black", stroke = 1), ncol = 1)) +
  theme(legend.text = element_text(size = 15))

PlotPDF('LCHADCellType.umap', 6, 8)
p1
dev.off()

p2 <- DimPlot2(LCHAD,
               reduction = 'umap.cca',
               group.by = 'Cell_state',
               cols = mycol_40,
               label = FALSE,
               raster = FALSE,
               pt.size = 0.2,
               order = TRUE,
               label.size = 10) +
  labs(title = 'Cell States', x = 'UMAP1', y = 'UMAP2') +
  guides(color = guide_legend(override.aes = list(size = 6, shape = 21, fill = mycol_40[1:32], color = "black", stroke = 1), ncol = 2)) +
  theme(legend.text = element_text(size = 15))

PlotPDF('LCHADCellStates.umap', 6, 8)
p2
dev.off()

## Nichenet Ligand/Receptor Analysis
library(nichenetr)

organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

Idents(LCHAD) <- LCHAD$group2

##Run Nichenet Wrapper on CM-->Myeloid
nichenet_outputCMtoMyeloid <- nichenet_seuratobj_aggregate(
  seurat_obj = LCHAD,
  sender = "CM",
  receiver = c("Myeloid"),
  condition_colname = "group2",
  condition_oi = "Peds_LCHAD",
  condition_reference = "Peds_Ctrl",
  expression_pct = 0.1,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks, top_n_ligands = 60, geneset = "up"
)

##Filter for DEGs
FinalCMbyGroupMarkers <- read_csv("FinalCMbyGroupMarkers.csv")
FinalCMPedLCHAD <- FinalCMbyGroupMarkers %>% filter(cluster == "Peds_LCHAD")
FinalLCHADMLLymphMarkers <- read_csv("FinalLCHADMLLymphMarkers.csv")
FinalMo-HLAhi <- FinalLCHADMLLymphMarkers %>% filter(cluster == "Mo-HLAhi")

nichenet_outputCMtoMyeloidLignandreceptors <- nichenet_outputCMtoMyeloid$ligand_receptor_df
filtered_listLigandtoreceptor <- nichenet_outputCMtoMyeloidLignandreceptors %>%
  filter(ligand %in% FinalCMPedLCHAD$gene & receptor %in% FinalMo-HLAhi$gene)

write.csv(filtered_listLigandtoreceptor, "CMtoMyeloidLigandtoreceptor.csv")

##Run Nichenet Wrapper on CM-->FB
nichenet_outputCMtoFB <- nichenet_seuratobj_aggregate(
  seurat_obj = LCHAD,
  sender = "CM",
  receiver = c("Fibroblast"),
  condition_colname = "group2",
  condition_oi = "Peds_LCHAD",
  condition_reference = "Peds_Ctrl",
  expression_pct = 0.1,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks, top_n_ligands = 60, geneset = "up"
)

##Filter for DEGs
FinalCMbyGroupMarkers <- read_csv("FinalCMbyGroupMarkers.csv")
FinalCMPedLCHAD <- FinalCMbyGroupMarkers %>% filter(cluster == "Peds_LCHAD")
FinalFBbyGroupMarkers <- read_csv("FinalFBbyGroupMarkers.csv")
FinalFBPedLCHAD <- FinalFBbyGroupMarkers %>% filter(cluster == "Peds_LCHAD")

nichenet_outputCMtoMyeloidLignandreceptors <- nichenet_outputCMtoFB$ligand_receptor_df
filtered_listLigandtoreceptor <- nichenet_outputCMtoMyeloidLignandreceptors %>%
  filter(ligand %in% FinalCMPedLCHAD$gene & receptor %in% FinalFBPedLCHAD$gene)

write.csv(filtered_listLigandtoreceptor, "CMtoFBLigandtoreceptor.csv")

####