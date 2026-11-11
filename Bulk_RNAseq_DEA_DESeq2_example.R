# Set global variables:
setwd("")
inputDir <- ""
outputDir <- ""
indexDir <- ""
treatment <- ""
Seq_Exp <- ""

# load packages
library(apeglm)
library(DESeq2)
library(pheatmap)
library(ggplot2)



# import data
countdata <- read.table(paste(inputDir, "", sep = ""), header = T)
gene_id <- read.csv(paste(indexDir, "", sep = ""), header = F)
gene_name <- read.csv(paste(indexDir, "", sep = ""), header = F)
gene_type <- read.csv(paste(indexDir, "", sep = ""), header = F)

## Data preparation
#---------------------------------------------------------------------
# how many lines 
nrow(gene_id)
nrow(gene_name)
nrow(gene_type)

# concatinate gene_id, gene_name, gene_type and make unique
gene_info <- data.frame (gene_id, gene_name, gene_type)
gene_info <- unique(gene_info)
nrow(gene_info)

# Delete and rename columns
countdata <- data.frame(countdata[,c(1,7:ncol(countdata))])
Newnames <- c("GeneID", gsub(".*\\.(.*?)_Aligned.*", "\\1", colnames(countdata[,c(2:ncol(countdata))])))
names(countdata) <- Newnames


# Rename gene info and make it unique
info.names <- c('GeneID', 'Gene_Name', 'Gene_Type')
names(gene_info) <- info.names

#sort both files before merge together
gene_info <- gene_info[order(gene_info$GeneID),]
countdata <- countdata[order(countdata$GeneID),]
# check if identical gene_id
all(gene_info$GeneID == countdata$GeneID)

# merge countdata and gene_info
countdata_transcript <- merge(gene_info,countdata, by="GeneID")

# gene_name as rownames
rownames(countdata_transcript) <- make.unique(as.character(countdata_transcript$Gene_Name))

# delete 3 columns
countdata_transcript <- countdata_transcript[-1]
countdata_transcript <- countdata_transcript[-1]
countdata_transcript <- countdata_transcript[-1]

countdata <- countdata_transcript[,c(5:8,1:4)]
############################################################################
# *Metadata
(condition <- factor(c("Ctrl", "Ctrl", "Ctrl", "Ctrl", rep(treatment, 4))))
(coldata <- data.frame(row.names=colnames(countdata),sample=colnames(countdata), condition))
#------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
## DESeq2 analysis
# DESeq object
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)

# Set factor and run DESeq
dds$condition <- factor(dds$condition, levels = c('Ctrl', treatment))
dds <- DESeq(dds, parallel = T)

# DESeq results object
res <- results(dds)
summary(res)
res

# Regularized log transformation for visualization
# Normalization
norm <- counts(dds, normalized=TRUE)
norm_dds <- norm

#rlogTransformation transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size
rld <- rlogTransformation(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef=paste("condition_",treatment,"_vs_Ctrl", sep = ""), type="apeglm")

# Adjusted p-value, up and down-regulated genes
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange > 0)
table(res$padj < 0.05 & res$log2FoldChange < 0)

# Print to table expressed, up and downregulated genes
up_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange > 0)
down_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange < 0)
diff_regulated <- subset(res, res$padj < 0.05)
all_genes <- res

write.table(up_regulated, file=paste(outputDir, Seq_Exp,"_upregulated_genes_LFC0.txt",sep = ""), sep = "\t", quote = FALSE)
write.table(down_regulated, file=paste(outputDir, Seq_Exp,"_downregulated_genes_LFC0.txt",sep = ""), sep = "\t", quote = FALSE)
write.table(diff_regulated, file=paste(outputDir, Seq_Exp, "_diffregulated_genes_LFC0.txt",sep = ""), sep = "\t", quote = FALSE)
write.table(all_genes, file=paste(outputDir, Seq_Exp, "_allregulated_genes_LFC0.txt",sep = ""), sep = "\t", quote = FALSE)

## Expressed genes
expressed_genes <- subset(res, res$baseMean > 10)
write.table(expressed_genes, file=paste(outputDir, Seq_Exp, "_expressed_genes.txt",sep = ""), sep = "\t", quote = FALSE)

# Visualization
# Principal Component Analysis - quality control
pdf(paste(outputDir, Seq_Exp, "_gene_PCA.pdf", sep=""))
plotPCA(rld) + theme_classic() + geom_label(aes(label = name), position = position_nudge(y = 1))
dev.off()


## Checking the validity of the cells by plotting expression of genes on a heatmap
## Check if the gene is present in the R-log object
pluripotency <- c("KLF4", "MYC", "NANOG", "LIN28A")
forebrain <- c("FOXG1", "OTX2", "PAX6", "EMX2", "NKX2-1", "DLX2", "LHX2")
mid_hindbrain <- c("OTX1", "EN1", "LMX1A", "CORIN", "FOXA2", "PAX8", "HOXA2")
neural <- c("SOX2", "TUBB", "DCX", "TBR1", "EOMES", "CUX2", "OLIG2", "GFAP", "EMX1", "NEUROD1")

gene_rld_df <- assay(rld)
biomarkers_count_norm <- gene_rld_df[which(rownames(gene_rld_df) %in% pluripotency | rownames(gene_rld_df) %in% forebrain | rownames(gene_rld_df) %in% mid_hindbrain | rownames(gene_rld_df) %in% neural),]
biomarkers_count_norm <- biomarkers_count_norm[c(pluripotency,
                                                 forebrain,
                                                 mid_hindbrain,
                                                 neural),]
coldata <- coldata[order(coldata$condition),]
annotation_row <- data.frame(gene_name = c(pluripotency,
                                           forebrain,
                                           mid_hindbrain,
                                           neural),
                             type = c(rep("pluripotency", length(pluripotency)),
                                      rep("forebrain", length(forebrain)),
                                      rep("mid_hindbrain", length(mid_hindbrain)),
                                      rep("neural", length(neural))))

rownames(annotation_row) <- annotation_row$gene_name
pdf(paste(outputDir, Seq_Exp, "_marker_genes_Heatmap.pdf", sep = ""))
pheatmap(biomarkers_count_norm[,coldata$sample], cluster_rows = F, cluster_cols = F, annotation_col = coldata[,"condition",drop=F], gaps_row = c(4,11,18), annotation_row = annotation_row[, "type", drop=F])
dev.off()
