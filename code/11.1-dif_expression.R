# Redirect all output to log file
con <- file(snakemake@log[[1]], "a+")
sink(con, append = TRUE, type="message")
sink(con, append = TRUE)

#=============================================================================#
# 0. Install packages if needed                                               #
#=============================================================================#
list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org');

to_install <- c("DESeq2", "ggplot2", "stringr", "apeglm", "pheatmap")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package);
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#
library("DESeq2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("stringr")
library(futile.logger)

#=============================================================================#
# 0. Logging                                                                  #
#=============================================================================#
# Create a logger that will be saved to a file
flog.logger("saturation", TRACE, appender=appender.file(snakemake@log[[1]]))

flog.info("R script to run a saturation analysis")

#=============================================================================#
# 1. Load input variables from Snakemake                                      #
#=============================================================================#
flog.info("Definining input variables")
sampleFiles <- snakemake@input[["counts"]]
sub_dir <- dirname(snakemake@output[["dif_expr"]][1])
plot_dir <- dirname(snakemake@output[["pca"]][1])
pcaplot <- snakemake@output[["pca"]]
distplot <- snakemake@output[["dist"]]
heatplot <- snakemake@output[["heatmap"]][1]

cw <- 18

# Make sure that output directories exist
if (!file.exists(sub_dir)){
  dir.create(sub_dir)
}

if (!file.exists(plot_dir)){
  dir.create(plot_dir)
}

#=============================================================================#
# 2. Read counts and create metadata                                          #
#=============================================================================#
flog.info("Reading input files")

# Assign a condition to each sample
sampleCondition <- c()
for (file in sampleFiles){
  if ((grepl("01", file, fixed = TRUE) == 1) | (grepl("02", file, fixed = TRUE) == 1)){
    sampleCondition <- c(sampleCondition, "mucoid")
  }
  else{sampleCondition <- c(sampleCondition, "inhibitor")}
}

# Create table with metadata
sampleTable <- data.frame(sampleName = substring(sampleFiles, first = 10, last = 11),
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

# Read counts from files
filename <- sampleFiles[1]
counts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
for (filename in sampleFiles[2:length(sampleFiles)]){
  cts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
  counts <- cbind(counts, cts)
}

# Add metadata (three conditions, mucoid/inhibitor, fructose/sucrose and date)
metadat <- snakemake@input[["meta"]]
coldata <- read.csv(metadat, sep = "\t", row.names = 1)
coldata$condition2 <- factor(coldata$condition2)
coldata$condition1 <- factor(coldata$condition1)
coldata$condition3 <- factor(coldata$condition3)

# Convert the data to DESeq input format
dds_pcr <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ condition3 + condition1 + condition2)

#=============================================================================#
# 3. Run a PCA for the complete dataset                                       #
#=============================================================================#
flog.info("Running PCA")
# PCA calculations
vsd <- vst(dds_pcr, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition1", "condition2"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plots and formatting
p <- ggplot(pcaData, aes(PC1, PC2, color=condition2, shape=condition1)) +
  geom_point(size=4, alpha=1, stroke=0) +
  guides(color = guide_legend(title = "Substrate"), shape = guide_legend(title = "Morphology")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed();
p <- p + stat_ellipse(geom="polygon", aes(fill = pcaData$condition2), 
  alpha = 0.2, show.legend = FALSE, level = 0.95) +
  theme_minimal() + theme(panel.grid = element_blank(), 
  panel.border = element_rect(fill= "transparent"));

# Saving the plots to files
ggsave(pcaplot[1], width = 9, height = 6);
ggsave(pcaplot[2], width = 9, height = 6);

#=============================================================================#
# 3. Create a global heatmap                                                  #
#=============================================================================#
flog.info("Creating a global heatmap")
dds_htmp <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ condition1 + condition2)
dds_htmp <- collapseReplicates(dds_htmp, as.factor(str_sub(colnames(dds_htmp), 1, 5)))
global_dds <- DESeq(dds_htmp)

res_htmp <- results(global_dds)
summary(res_htmp);
res_htmp_filter <- subset(res_htmp, padj < .1);
res_htmp_filter <- subset(res_htmp_filter, baseMean >= 100);

vsd_htmp <- vst(dds_htmp, blind=FALSE)

htmp_select <- res_htmp_filter[order(abs(res_htmp_filter$log2FoldChange),
                                    decreasing = TRUE),]
htmp_select <- rownames(htmp_select)[1:30]
select <- rownames(dds_htmp) %in% htmp_select

df <- as.data.frame(colData(global_dds)[,c("condition1", "condition2")])
colnames(df) <- c("Phenotype", "Carbon source")
rownames(vsd_htmp) <- str_sub(rownames(vsd_htmp), 4, -1)
rownames(df) <- sub("_", "", rownames(df))
colnames(vsd_htmp) <- sub("_", "", colnames(vsd_htmp))

png(file = heatplot, width = 600, height = 400);
pheatmap(assay(vsd_htmp)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=40)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][2], width = 600, height = 400);
pheatmap(assay(vsd_htmp)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=40)
invisible(dev.off())

sampleDists <- dist(t(assay(vsd_htmp)))

#=============================================================================#
# 4. Create a distance plot                                                   #
#=============================================================================#
flog.info("Creating a distance plot")
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "BuPu")) )(255)

png(file = distplot[1], width = 600, height = 400);
postscript(file = distplot[2], width = 600, height = 400);
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
invisible(dev.off())

#=============================================================================#
# 5. Differential expression analyses                                         #
#=============================================================================#
# Run a differential expression analysis for Smucoid vs Sinhibitor
flog.info("Running DESeq2")
savename <- "Smucoid_vs_Sinhibitor"
flog.info(paste("Comparison", savename, sep = " "))

# Get data for the differential expression analysis
coldata2 <- coldata[coldata$condition2 == "S",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition1)

# Run the differential expression analysis
dds <- DESeq(dds)


# Get results with and without lfc threshold filters
res <- lfcShrink(dds, coef="condition1_mucoid_vs_inhibitor", type="apeglm")
res_subset <- results(dds, lfcThreshold=1);
summary(res);

res_subset <- subset(res_subset, padj < .1);
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100);

# Create a heatmap for this comparison
res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),]
res_select <- rownames(res_select)[1:20]
select <- rownames(dds) %in% res_select

df <- as.data.frame(colData(dds)[,c("condition1", "condition3")])
colnames(df) <- c("Phenotype", "Batch")
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = "_S_")]
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1)

png(file = snakemake@output[["heatmap"]][3], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][4], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

# Save results to csv files
write.csv(res_filter, file = snakemake@output[["dif_expr"]][1], quote=FALSE);
write.csv(res_subset, file = snakemake@output[["dif_expr"]][2], quote=FALSE);

# Generate MA plots (png and postscript)
png(file = snakemake@output[["plots"]][1], width = 600, height = 400);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())
postscript(file = snakemake@output[["plots"]][2], width = 900, height = 600);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())

# Run a differential expression analysis for Fmucoid vs Finhibitor
savename <- "Fmucoid_vs_Finhibitor"
flog.info(paste("Comparison", savename, sep = " "))
coldata2 <- coldata[coldata$condition2 == "F",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition1)

dds <- DESeq(dds)

res <- lfcShrink(dds, coef="condition1_mucoid_vs_inhibitor", type="apeglm")
res_subset <- results(dds, lfcThreshold=1);
summary(res);

res_subset <- subset(res_subset, padj < .1);
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100);

res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),]
res_select <- rownames(res_select)[1:20]
select <- rownames(dds) %in% res_select

df <- as.data.frame(colData(dds)[,c("condition1", "condition3")])
colnames(df) <- c("Phenotype", "Batch")
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = "_F_")]
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1)

png(file = snakemake@output[["heatmap"]][5], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][6], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

write.csv(res_filter, file = snakemake@output[["dif_expr"]][3], quote=FALSE);
write.csv(res_subset, file = snakemake@output[["dif_expr"]][4], quote=FALSE);

postscript(file=snakemake@output[["plots"]][4], width=900, height=600);
plotMA(res, ylim=c(-3,3), colSig = "#c00000");
invisible(dev.off())
png(file=snakemake@output[["plots"]][3], width=600, height=400);
plotMA(res, ylim=c(-3,3), colSig = "#c00000");
invisible(dev.off())


# Run a differential expression analysis for Smucoid vs Fmucoid
savename <- "Smucoid_vs_Fmucoid"
flog.info(paste("Comparison", savename, sep = " "))
coldata2 <- coldata[coldata$condition1 == "mucoid",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition2)

dds <- DESeq(dds)

res <- lfcShrink(dds, coef="condition2_S_vs_F", type="apeglm")
res_subset <- results(dds, lfcThreshold=1);
summary(res);

res_subset <- subset(res_subset, padj < .1);
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100);

res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),]
res_select <- rownames(res_select)[1:20]
select <- rownames(dds) %in% res_select

df <- as.data.frame(colData(dds)[,c("condition2", "condition3")])
colnames(df) <- c("Phenotype", "Batch")
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = "I01|I02")]
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1)

png(file = snakemake@output[["heatmap"]][7], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][8], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

write.csv(res_filter, file = snakemake@output[["dif_expr"]][5], quote=FALSE);
write.csv(res_subset, file = snakemake@output[["dif_expr"]][6], quote=FALSE);

png(file = snakemake@output[["plots"]][5], width = 600, height = 400);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())
postscript(file = snakemake@output[["plots"]][6], width = 900, height = 600);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())

# Run a differential expression analysis for Sinhibitor vs Finhibitor
savename <- "Sinhibitor_vs_Finhibitor"
flog.info(paste("Comparison", savename, sep = " "))
coldata2 <- coldata[coldata$condition1 == "inhibitor",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition2)

dds <- DESeq(dds)

res <- lfcShrink(dds, coef="condition2_S_vs_F", type="apeglm")
res_subset <- results(dds, lfcThreshold=1);
summary(res);

res_subset <- subset(res_subset, padj < .1);
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100);

res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),]
res_select <- rownames(res_select)[1:20]
select <- rownames(dds) %in% res_select

df <- as.data.frame(colData(dds)[,c("condition2", "condition3")])
colnames(df) <- c("Phenotype", "Batch")
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = "I09|I10")]
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1)

png(file = snakemake@output[["heatmap"]][9], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][10], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

write.csv(res_filter, file = snakemake@output[["dif_expr"]][7], quote=FALSE);
write.csv(res_subset, file = snakemake@output[["dif_expr"]][8], quote=FALSE);

png(file = snakemake@output[["plots"]][7], width = 600, height = 400);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())
postscript(file = snakemake@output[["plots"]][8], width = 900, height = 600);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())