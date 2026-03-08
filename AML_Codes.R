# =========================================================
# AML scRNA-seq analysis
# Personal paths removed / minimally cleaned for GitHub
# =========================================================

# ---------------------------
# Libraries
# ---------------------------
library(GSVA)
library(ggrepel)
library(monocle3)
library(Seurat)
library(Nebulosa)
library(CellChat)
library(ggpubr)
library(scDblFinder)
library(scRNAseq)
library(dplyr)
library(SingleR)
library(ggplot2)
library(dittoSeq)
library(escape)
library(scales)
library(RColorBrewer)
library(viridis)
library(harmony)
library(scRepertoire)

# Optional packages used later
library(GSEABase)
library(Polychrome)
library(pals)
library(wesanderson)

# ---------------------------
# Paths
# ---------------------------
# Replace these with your local/project folders
data_dir <- "data"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Read 10X data: Controls
# ---------------------------
HC1.data <- Read10X(data.dir = file.path(data_dir, "C1"))
HC2.data <- Read10X(data.dir = file.path(data_dir, "C2"))
HC3.data <- Read10X(data.dir = file.path(data_dir, "C3"))
HC4.data <- Read10X(data.dir = file.path(data_dir, "C4"))
HC5.data <- Read10X(data.dir = file.path(data_dir, "C5"))
HC6.data <- Read10X(data.dir = file.path(data_dir, "C6"))
HC7.data <- Read10X(data.dir = file.path(data_dir, "C7"))

HC1 <- CreateSeuratObject(counts = HC1.data$`Gene Expression`, project = "HC1")
HC2 <- CreateSeuratObject(counts = HC2.data$`Gene Expression`, project = "HC2")
HC3 <- CreateSeuratObject(counts = HC3.data$`Gene Expression`, project = "HC3")
HC4 <- CreateSeuratObject(counts = HC4.data, project = "HC4")
HC5 <- CreateSeuratObject(counts = HC5.data$`Gene Expression`, project = "HC5")
HC6 <- CreateSeuratObject(counts = HC6.data$`Gene Expression`, project = "HC6")
HC7 <- CreateSeuratObject(counts = HC7.data$`Gene Expression`, project = "HC7")

HC1$Subset <- "Control"
HC2$Subset <- "Control"
HC3$Subset <- "Control"
HC4$Subset <- "Control"
HC5$Subset <- "Control"
HC6$Subset <- "Control"
HC7$Subset <- "Control"

# ---------------------------
# Read 10X data: AML
# ---------------------------
P1.data <- Read10X(data.dir = file.path(data_dir, "P1"))
P2.data <- Read10X(data.dir = file.path(data_dir, "P2"))
P3.data <- Read10X(data.dir = file.path(data_dir, "P3"))
P4.data <- Read10X(data.dir = file.path(data_dir, "P4"))
P5.data <- Read10X(data.dir = file.path(data_dir, "P5"))
P6.data <- Read10X(data.dir = file.path(data_dir, "P6"))
P7.data <- Read10X(data.dir = file.path(data_dir, "P7"))
P8.data <- Read10X(data.dir = file.path(data_dir, "P8"))
P9.data <- Read10X(data.dir = file.path(data_dir, "P9"))
P10.data <- Read10X(data.dir = file.path(data_dir, "P10"))

P1 <- CreateSeuratObject(counts = P1.data, project = "P1")
P2 <- CreateSeuratObject(counts = P2.data, project = "P2")
P3 <- CreateSeuratObject(counts = P3.data, project = "P3")
P4 <- CreateSeuratObject(counts = P4.data, project = "P4")
P5 <- CreateSeuratObject(counts = P5.data, project = "P5")
P6 <- CreateSeuratObject(counts = P6.data, project = "P6")
P7 <- CreateSeuratObject(counts = P7.data, project = "P7")
P8 <- CreateSeuratObject(counts = P8.data, project = "P8")
P9 <- CreateSeuratObject(counts = P9.data, project = "P9")
P10 <- CreateSeuratObject(counts = P10.data, project = "P10")

P1$Subset <- "AML"
P2$Subset <- "AML"
P3$Subset <- "AML"
P4$Subset <- "AML"
P5$Subset <- "AML"
P6$Subset <- "AML"
P7$Subset <- "AML"
P8$Subset <- "AML"
P9$Subset <- "AML"
P10$Subset <- "AML"

# ---------------------------
# Merge objects
# ---------------------------
JSS.BIG <- merge(
  HC1,
  y = c(HC2, HC3, HC4, HC5, HC6, HC7, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10),
  add.cell.ids = c(
    "HC1", "HC2", "HC3", "HC4", "HC5", "HC6", "HC7",
    "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"
  ),
  project = "AML"
)

# ---------------------------
# QC
# ---------------------------
JSS.BIG[["percent.mt"]] <- PercentageFeatureSet(JSS.BIG, pattern = "^MT-")

VlnPlot(
  JSS.BIG,
  features = c("nFeature_RNA", "percent.mt"),
  group.by = "Subset",
  ncol = 3,
  pt.size = 0
)
ggsave(file.path(output_dir, "Vlnplot_QC.pdf"), width = 15, height = 10)

# ---------------------------
# Integration
# ---------------------------
JSS.BIG <- JoinLayers(JSS.BIG)
JSS.BIG <- subset(JSS.BIG, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

JSS_list <- SplitObject(JSS.BIG, split.by = "orig.ident")
JSS_list <- lapply(X = JSS_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = JSS_list)
JSS.anchors <- FindIntegrationAnchors(object.list = JSS_list, anchor.features = features)
JSS.combined <- IntegrateData(anchorset = JSS.anchors)

DefaultAssay(JSS.combined) <- "integrated"
JSS.combined <- ScaleData(JSS.combined, verbose = FALSE)
JSS.combined <- RunPCA(JSS.combined, npcs = 30, verbose = FALSE)
JSS.combined <- RunUMAP(JSS.combined, reduction = "pca", dims = 1:30)
JSS.combined <- RunTSNE(JSS.combined, reduction = "pca", dims = 1:30)
JSS.combined <- FindNeighbors(JSS.combined, reduction = "pca", dims = 1:30)
JSS.combined <- FindClusters(JSS.combined, resolution = 1.2)

DefaultAssay(JSS.combined) <- "RNA"
JSS.combined <- ScaleData(JSS.combined, verbose = FALSE)
saveRDS(JSS.combined, file.path(output_dir, "AML.combined.rds"))

JSS.combined <- JoinLayers(JSS.combined)

# ---------------------------
# Basic plots
# ---------------------------
DimPlot(JSS.combined)

# ---------------------------
# Gene set preparation
# ---------------------------
sign <- read.delim(file.path(data_dir, "immune_sign_Immune_final.txt"))

full <- as.list(sign)
unique_names <- names(full)
list <- list()

for (i in seq_along(unique_names)) {
  tmp <- full[[i]]
  tmp <- tmp[tmp != ""]
  tmp <- unique(toupper(tmp))
  tmp <- GSEABase::GeneSet(tmp, setName = paste(unique_names[i]))
  list[[i]] <- tmp
}
list <- GSEABase::GeneSetCollection(list)

ES.immune <- escape.matrix(
  JSS.combined,
  gene.sets = list,
  min.size = 3,
  method = "ssGSEA"
)
saveRDS(ES.immune, file.path(output_dir, "ES.immune.rds"))

# ---------------------------
# Visualization
# ---------------------------
DimPlot(JSS.combined, group.by = "Subset", reduction = "tsne", label = FALSE) + NoLegend()
ggsave(file.path(output_dir, "Dimplot_Ctrl_AML.svg"), width = 5, height = 5)

DimPlot(JSS.combined, reduction = "tsne", label = TRUE) + NoLegend()
ggsave(file.path(output_dir, "Dimplot.svg"), width = 5, height = 5)

FeaturePlot(JSS.combined, features = c("CD34"), label = TRUE, raster = FALSE, order = TRUE, reduction = "tsne")
ggsave(file.path(output_dir, "Featureplot_CD34.png"), width = 5, height = 5)

FeaturePlot(JSS.combined, features = c("CD34"), raster = FALSE, label = TRUE, order = TRUE, reduction = "tsne")
ggsave(file.path(output_dir, "Featureplot_CD34.svg"), width = 5, height = 5)

DotPlot(JSS.combined, features = c("CD34"))
ggsave(file.path(output_dir, "Dotplot_CD34.svg"), width = 4, height = 6)

FeaturePlot(JSS.combined, features = c("IKBKB"), raster = FALSE, label = TRUE, order = TRUE, reduction = "tsne")
ggsave(file.path(output_dir, "Featureplot_IKBKB.png"), width = 5, height = 5)

# ---------------------------
# IKBKB / NR4A1 subset analysis
# ---------------------------
JSS.combined$Ikbkb <- JSS.combined@assays$RNA$data["IKBKB", ]
JSS.combined$IkbkbSep <- ifelse(JSS.combined$Ikbkb > 0, "Pos", "Neg")
JSS.combined.Ikbkb <- subset(JSS.combined, IkbkbSep %in% c("Pos"))

JSS.combined.Ikbkb$Nr4a1Sep <- ifelse(JSS.combined.Ikbkb$Nr4a1 > 0, "Pos", "Neg")
JSS.combined.IkbkbNr4a1 <- subset(JSS.combined.Ikbkb, Nr4a1Sep %in% c("Pos"))

DimPlot(JSS.combined.IkbkbNr4a1, label = TRUE)
ggsave(file.path(output_dir, "Dimplot_DP.png"), width = 7, height = 5)

FeaturePlot(JSS.combined.IkbkbNr4a1, features = c("NR4A1", "IKBKB"))
ggsave(file.path(output_dir, "Featureplot_DP_IKBKB_NR4A1.png"), width = 7, height = 5)

FeatureScatter(
  subset(JSS.combined.IkbkbNr4a1, Subset %in% c("HC") & seurat_clusters %in% c(0)),
  feature1 = c("NR4A1"),
  feature2 = c("IKBKB"),
  jitter = TRUE,
  log = TRUE
)
ggsave(file.path(output_dir, "FeatureScatter_DP_AML_C26_IKBKB_NR4A1.png"), width = 4, height = 4)

DotPlot(JSS.combined, features = c("NR4A1"))
DotPlot(JSS.combined.IkbkbNr4a1, features = c("NR4A1", "IKBKB"))

FeatureScatter(
  JSS.combined.IkbkbNr4a1,
  feature1 = "NR4A1",
  feature2 = c("IKBKB"),
  log = TRUE,
  jitter = TRUE,
  span = TRUE
)
ggsave(file.path(output_dir, "Featureplot_BothPositive_IKBKB_NR4A1.png"), width = 7, height = 5)

FeatureScatter(
  subset(JSS.combined.IkbkbNr4a1, seurat_clusters %in% c(0, 4, 12, 13, 24, 26, 29, 30)),
  feature1 = "NR4A1",
  feature2 = c("IKBKB"),
  log = TRUE,
  jitter = TRUE,
  span = TRUE
)
ggsave(file.path(output_dir, "Featureplot_BothPositive_AMLclusters_IKBKB_NR4A1.png"), width = 7, height = 5)

table(JSS.combined.IkbkbNr4a1$seurat_clusters, JSS.combined.IkbkbNr4a1$Subset)