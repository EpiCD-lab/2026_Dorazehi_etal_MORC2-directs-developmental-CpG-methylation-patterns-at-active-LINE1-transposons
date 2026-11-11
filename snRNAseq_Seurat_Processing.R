## ==== The following code takes the output from cellRanger and produces the plots shown in the paper ==== ##
## Code for running cellRanger:
## << cellranger count --include-introns --id sample_id --sample sample_id  --transcriptome path/to/cellRanger/index --fastqs /path/to/fastq/dir >>


##==== Rcode of an example analysis ====##

## Load libraries
library(SingleCellExperiment)
library(harmony)
library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(Seurat)
library(stringr)
library(RColorBrewer)


## load global variables
setwd("Projects/MORC2/")
input_prefix <- "snRNAseq_analysis/"
input_suffix <- "/outs/filtered_feature_bc_matrix/"
plots_outDir <- "Results/snRNAseq_out/plots/"
data_outDir <- "Results/snRNAseq_out/data/"
rds_outDir <- "Results/snRNAseq_out/rds/"
experiment <- "h6_WT_vs_MORC2_CRISPRi"

####################################################################################################
## SAMPLES:
## snRNAseq_d15org_H6_LacZ ##
## snRNAseq_d15org_H6_MORC2i ##
####################################################################################################

## Load samples from cellRanger output:
samp.list <- list.files(path = paste0(input_prefix), full.names = F, all.files = F, include.dirs = T)
samp.list <- samp.list[grep("snRNAseq_d15org_H6", samp.list)]

for (sample in samp.list) {
  ## Loading data and extracting raw counts --> use filtered_feature_bc_matrix
  control.counts <- Read10X(data.dir = paste(input_prefix, sample, input_suffix, sep = ""))
  control.counts <- CreateSeuratObject(counts = control.counts, project = sample, min.cells = 3, min.features = 200)
  control.counts
  
  ## QC ##
  # Percent of mitochondrial and ribosomal genes:
  control.counts[["percent.mt"]] <- PercentageFeatureSet(control.counts, pattern = "^MT-")
  control.counts[["percent.rb"]] <- PercentageFeatureSet(control.counts, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  control.counts$log10GenesPerUMI <- log10(control.counts$nFeature_RNA) / log10(control.counts$nCount_RNA)
  
  # Visualize QC metrics as a violin plot
  pdf(paste(plots_outDir, sample, "nFeat_RNA_nCount_RNA_preFiltering_QC.pdf", sep = ""))
  VlnPlot(control.counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  ## Filtering data according to gene expression, number of genes expressed and mitochrondrial genes expressed:
  upr_barcode_lmt <- mean(control.counts@meta.data$nFeature_RNA) + 2*sd(control.counts@meta.data$nCount_RNA)
  lwr_barcode_lmt <- mean(control.counts@meta.data$nFeature_RNA) - 2*sd(control.counts@meta.data$nFeature_RNA)
  
  control.counts <- subset(control.counts, subset = nFeature_RNA > 600 &
                             nFeature_RNA < upr_barcode_lmt &
                             percent.mt < 2)

  # Visualize QC as a violin plot post QC
  pdf(paste(plots_outDir, sample, "nFeat_RNA_nCount_RNA_postFiltering_QC.pdf", sep = ""))
  VlnPlot(control.counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  ## Normalizing data:
  control.counts.norm <- NormalizeData(control.counts, normalization.method = "LogNormalize", scale.factor = 10000)
  control.counts.norm <- NormalizeData(control.counts.norm)
  
  ## SAVE RDS files here and continue with them together
  rds.save.list <- list(control.counts, control.counts.norm)
  saveRDS(rds.save.list, paste(rds_outDir, sample, ".rds", sep = ""))
  print(paste0("Sample ", which(samp.list == sample), ": done"))
}

############################################################
############################################################
## Created RDS file containing counts + normalized_counts for Controls and Treatments
# Loading the RDS files:
rds_files <- list.files(path = rds_outDir, pattern = "\\.rds$", full.names = TRUE)
rds_files <- rds_files[grep("snRNAseq_d15org_H6", rds_files)] ## LacZ vs MORC2 CRISPRi

# Load each RDS file into a list
loaded_objects <- lapply(rds_files, readRDS)
names(loaded_objects) <- basename(rds_files)

## Adding metadata [NOTE!! Adjust accordingly to the data that you are using]:
for (sample in names(loaded_objects)) {
  loaded_objects[[sample]][[2]]@meta.data$condition <- sapply(strsplit(sample, split = "_"), function(x) ifelse(length(x) >= 5, x[2], NA))
  #loaded_objects[[sample]][[2]]@meta.data$replicate <- sapply(strsplit(sample, split = "_"), function(x) ifelse(length(x) >= 5, x[3], NA))
  #loaded_objects[[sample]][[2]]@meta.data$batch <- gsub(".*_(batch[0-9]).*", "\\1", sample)
}

MORC2.CRISPRi <- lapply(loaded_objects, function(x) x[[2]])
MORC2.CRISPRi <- Reduce(function(x, y) merge(x, y, add.cell.ids = c("Object1", "Object2")), MORC2.CRISPRi)


MORC2.CRISPRi <- FindVariableFeatures(MORC2.CRISPRi, selection.method = "vst", nfeatures = 2000)

## All samples post filtering QC
pdf(paste0(plots_outDir, experiment,  "_nFeat_RNA_nCount_RNA_postFiltering_QC.pdf", sep = ""))
VlnPlot(MORC2.CRISPRi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
dev.off()

## Scaling the data (this step is standard. Not sure if only a subset of cells or all?)
MORC2.CRISPRi <- ScaleData(MORC2.CRISPRi)
MORC2.CRISPRi <- RunPCA(TRIM28.CRISPRi, features = VariableFeatures(object = MORC2.CRISPRi))

## Use when needed
#set.seed(123)
#MORC2.CRISPRi.harmony <- RunHarmony(MORC2.CRISPRi,
#                                      c("batch", "replicate"),
#                                      plot_convergence = TRUE)


# Examine and visualize PCA results a few different ways
print(MORC2.CRISPRi[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(MORC2.CRISPRi, dims = 1:2, reduction = "umap")

DimPlot(MORC2.CRISPRi, reduction = "umap", split.by = "condition")
DimPlot(MORC2.CRISPRi, reduction = "umap", group.by = "condition")
DimPlot(MORC2.CRISPRi, reduction = "umap", split.by = "replicate")
DimPlot(MORC2.CRISPRi, reduction = "umap", split.by = "batch")

pdf(paste0(plots_outDir, experiment,  "_umap_perSample.pdf", sep = ""), width = 60)
DimPlot(MORC2.CRISPRi, reduction = "umap", split.by = "orig.ident")
dev.off()

DimHeatmap(object = MORC2.CRISPRi, reduction = "umap", cells = 500, dims = 1:10)
ElbowPlot(MORC2.CRISPRi, ndims = 50, reduction = "umap")

MORC2.CRISPRi <- FindNeighbors(MORC2.CRISPRi, dims = 1:30, reduction = "umap") ## alter as needed
MORC2.CRISPRi <- FindClusters(MORC2.CRISPRi, resolution = 0.1) ## alter as needed

## Establishing cell population identity:
pdf(paste0(plots_outDir, experiment,  "_umap.pdf", sep = ""))
DimPlot(MORC2.CRISPRi, reduction = "umap")
dev.off()

pdf(paste0(plots_outDir, experiment,  "_umap_conditionSplit.pdf", sep = ""))
DimPlot(MORC2.CRISPRi, reduction = "umap", split.by = "condition")
dev.off()

pdf(paste0(plots_outDir, experiment,  "_cluster_gene_markers.pdf", sep = ""))
DotPlot(MORC2.CRISPRi, c("GAD1", "DCX", "TBR1", "RBFOX3", "NCAM1",
                         "MAP2", "BCL11B", "EOMES", "PCNA", "SOX2",
                         "PAX6", "CUX2", "GLI3", "TJP1",
                         "POU5F1", "MYC", "KLF4"),
        cols = c("white", "#005A32")) + coord_flip()
dev.off()

FeaturePlot(MORC2.CRISPRi, "MORC2", order = T, max.cutoff = 2, split.by = "condition")

## Define percentage of cells in population
MORC2.CRISPRi.relab <- c("NPC", "NPC", "Others", "Early_neurons", "Early_neurons")


for (cellType in 1:length(MORC2.CRISPRi.relab)) {
  MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == cellType-1)] <- MORC2.CRISPRi.relab[cellType]
}

#MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == 0)] <- "NPC"
#MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == 1)] <- "NPC"
#MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == 2)] <- "Others"
#MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == 3)] <- "Early_neurons"
#MORC2.CRISPRi@meta.data$seurat_clusters_relabeled[which(MORC2.CRISPRi@meta.data$seurat_clusters == 4)] <- "Early_neurons"

DimPlot(MORC2.CRISPRi, reduction = "umap", group.by = "seurat_clusters_relabeled", split.by = "condition")

cluster.cell.identity <- data.frame(clusters_relabeled = MORC2.CRISPRi@meta.data$seurat_clusters_relabeled,
                                    orig.identity = MORC2.CRISPRi@meta.data$orig.ident,
                                    clusters = MORC2.CRISPRi@meta.data$seurat_clusters,
                                    condition = MORC2.CRISPRi@meta.data$condition)
                                    #experiment = MORC2.CRISPRi@meta.data$replicate, ## Use if needed
                                    #batch = MORC2.CRISPRi@meta.data$batch) ## Use if needed

data_summary <- cluster.cell.identity %>%
  group_by(clusters_relabeled, condition) %>%
  summarise(count = n()) %>%
  group_by(condition) %>%
  mutate(percentage = count / sum(count) * 100)

ggplot(data_summary, aes(x = condition, y = percentage, fill = clusters_relabeled)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

data_summary <- cluster.cell.identity %>%
  group_by(clusters, experiment) %>%
  summarise(count = n()) %>%
  group_by(experiment) %>%
  mutate(percentage = count / sum(count) * 100)

ggplot(data_summary, aes(x = experiment, y = percentage, fill = clusters)) +
  geom_bar(stat = 'identity', position = 'stack', color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()
