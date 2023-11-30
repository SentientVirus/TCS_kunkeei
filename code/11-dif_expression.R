#=============================================================================#
# 0. Install packages if needed                                               #
#=============================================================================#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org')

to_install <- c("DESeq2", "ggplot2", "stringr")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package)
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#
library("DESeq2")
library("ggplot2")
library("stringr")

#=============================================================================#
# 1. Load input variables from Snakemake                                      #
#=============================================================================#
sampleFiles <- snakemake@input[["counts"]]
sub_dir <- dirname(snakemake@output[["DE"]][1])
plot_dir <- dirname(snakemake@output[["PCA"]][1])
pcaplot <- snakemake@output[["PCA"]]

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
  coord_fixed()
p <- p + stat_ellipse(geom="polygon", aes(fill = pcaData$condition2), 
  alpha = 0.2, show.legend = FALSE, level = 0.95) +
  theme_minimal() + theme(panel.grid = element_blank(), 
  panel.border = element_rect(fill= "transparent"))

# Saving the plots to files
ggsave(pcaplot[1], width = 9, height = 6)
ggsave(pcaplot[2], width = 9, height = 6)

#=============================================================================#
# 4. Differential expression analyses                                         #
#=============================================================================#
# Get data for the differential expression analysis Smucoid vs Sinhibitor
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
res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

# Save results to csv files
savename <- "Smucoid_vs_Sinhibitor"
write.csv(res_filter, file = snakemake@output[["DE"]][1], quote=FALSE)
write.csv(res_subset, file = snakemake@output[["DE"]][2], quote=FALSE)

# Generate MA plots (png and postscript)
png(file = snakemake@output[["plots"]][1], width = 600, height = 400)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()
postscript(file = snakemake@output[["plots"]][2], width = 900, height = 600)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()

# Run a differential expression analysis for Fmucoid vs Finhibitor
coldata2 <- coldata[coldata$condition2 == "F",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition1)

dds <- DESeq(dds)

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Fmucoid_vs_Finhibitor"
write.csv(res_filter, file = snakemake@output[["DE"]][3], quote=FALSE)
write.csv(res_subset, file = snakemake@output[["DE"]][4], quote=FALSE)

postscript(file=psnakemake@output[["plots"]][4], width=900, height=600)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
png(file=snakemake@output[["plots"]][3], width=600, height=400)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()


# Run a differential expression analysis for Smucoid vs Fmucoid
coldata2 <- coldata[coldata$condition1 == "mucoid",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition2)

dds <- DESeq(dds)

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Smucoid_vs_Fmucoid"
write.csv(res_filter, file = snakemake@output[["DE"]][5], quote=FALSE)
write.csv(res_subset, file = psnakemake@output[["DE"]][6], quote=FALSE)

png(file = snakemake@output[["plots"]][5], width = 600, height = 400)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()
postscript(file = snakemake@output[["plots"]][6], width = 900, height = 600)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()

# Run a differential expression analysis for Sinhibitor vs Finhibitor
coldata2 <- coldata[coldata$condition1 == "inhibitor",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition2)

dds <- DESeq(dds)

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Sinhibitor_vs_Finhibitor"
write.csv(res_filter, file = snakemake@output[["DE"]][7], quote=FALSE)
write.csv(res_subset, file = snakemake@output[["DE"]][8], quote=FALSE)

png(file = snakemake@output[["plots"]][7], width = 600, height = 400)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()
postscript(file = snakemake@output[["plots"]][8], width = 900, height = 600)
plotMA(res, ylim = c(-3,3), colSig = "#c00000")
dev.off()
