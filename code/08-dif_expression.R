if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

to_install <- c("DESeq2", "ggplot2", "stringr")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package)
}

library("DESeq2")
library("ggplot2")
library("stringr")

meta_dir <- "featureCounts_reverse/countfiles"
directory <- "featureCounts_reverse/countfiles/filtered"
sub_dir <- "results/DE"
plot_dir <- "plots"

# Make sure that output directories exist
if (!file.exists(sub_dir)){
  dir.create(sub_dir)
}

if (!file.exists(plot_dir)){
  dir.create(plot_dir)
}

# Assign a condition to each sample
sampleFiles <- list.files(directory)[grepl("counts.tsv",list.files(directory))][-1]
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
filename <- sprintf("%s/%s", directory, sampleFiles[1])
counts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
for (file in sampleFiles[2:length(sampleFiles)]){
  filename <- sprintf("%s/%s", directory, file)
  cts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
  counts <- cbind(counts, cts)
}

# Add metadata (three conditions, mucoid/inhibitor, fructose/sucrose and date)
metadat <- grep("metadata.tsv", list.files(meta_dir), value = TRUE)
coldata <- read.csv(sprintf("%s/%s", meta_dir, metadat), sep = "\t", row.names = 1)
coldata$condition2 <- factor(coldata$condition2)
coldata$condition1 <- factor(coldata$condition1)
coldata$condition3 <- factor(coldata$condition3)

# Convert the data to DESeq input format
dds_pcr <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ condition3 + condition1 + condition2)

# Run a PCA for the complete dataset
pdf(file = paste(plot_dir,"pcaplot.pdf", sep = "/"), width = 900, height = 600)
vsd <- vst(dds_pcr, blind=FALSE)
plotPCA(vsd, intgroup=c("condition1", "condition2"))
pcaData <- plotPCA(vsd, intgroup=c("condition1", "condition2"), returnData=TRUE)
substrate <- pcaData$condition2
morphology <- pcaData$condition1
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition2, shape=condition1)) +
  geom_point(size=4, alpha=1, stroke=0) +
  guides(color = guide_legend(title = "Substrate"), shape = guide_legend(title = "Morphology")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
p <- p + stat_ellipse(geom="polygon", aes(fill = pcaData$condition2), 
                      
                      alpha = 0.2, 
                      
                      show.legend = FALSE, 
                      
                      level = 0.95) +
  
  theme_minimal() +
  
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent"))

p
ggsave(paste(plot_dir,"pcaplot.png", sep = "/"), width = 9, height = 6)

# Run a differential expression analysis for Smucoid vs Sinhibitor
coldata2 <- coldata[coldata$condition2 == "S",]
coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
coldata2$condition3 = droplevels(coldata2$condition3)

dds <- DESeqDataSetFromMatrix(countData = counts2,
                              colData = coldata2,
                              design = ~ condition3 + condition1)

dds <- DESeq(dds)

#keep <- rowSums(counts(dds) >= 1) >= 3
#dds <- dds[keep,]

# Get results with and without lfc threshold filters
res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

# Save results to csv files
savename <- "Smucoid_vs_Sinhibitor"
write.csv(res_filter, file = paste(sub_dir, paste(savename, "csv", sep = "."), sep ="/"), quote=FALSE)
write.csv(res_subset, file = paste(sub_dir, paste(savename, "lfc1.csv", sep = "_"), sep ="/"), quote=FALSE)

# Generate MA plots (png and postscript)
postscript(file=paste(plot_dir, paste(savename, "ps", sep = "."), sep ="/"), width=900, height=600)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
png(file=paste(plot_dir, paste(savename, "png", sep = "."), sep ="/"), width=600, height=400)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
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

#keep <- rowSums(counts(dds) >= 1) >= 3
#dds <- dds[keep,]

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Fmucoid_vs_Finhibitor"
write.csv(res_filter, file = paste(sub_dir, paste(savename, "csv", sep = "."), sep ="/"), quote=FALSE)
write.csv(res_subset, file = paste(sub_dir, paste(savename, "lfc1.csv", sep = "_"), sep ="/"), quote=FALSE)

postscript(file=paste(plot_dir, paste(savename, "ps", sep = "."), sep ="/"), width=900, height=600)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
png(file=paste(plot_dir, paste(savename, "png", sep = "."), sep ="/"), width=600, height=400)
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

#keep <- rowSums(counts(dds) >= 1) >= 3
#dds <- dds[keep,]

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Smucoid_vs_Fmucoid"
write.csv(res_filter, file = paste(sub_dir, paste(savename, "csv", sep = "."), sep ="/"), quote=FALSE)
write.csv(res_subset, file = paste(sub_dir, paste(savename, "lfc1.csv", sep = "_"), sep ="/"), quote=FALSE)


postscript(file=paste(plot_dir, paste(savename, "ps", sep = "."), sep ="/"), width=900, height=600)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
png(file=paste(plot_dir, paste(savename, "png", sep = "."), sep ="/"), width=600, height=400)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
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

#keep <- rowSums(counts(dds) >= 1) >= 3
#dds <- dds[keep,]

res <- results(dds)
res_subset <- results(dds, lfcThreshold=1)
res
res_filter <- subset(res, padj < .1)
res_subset <- subset(res_subset, padj < .1)

savename <- "Sinhibitor_vs_Finhibitor"
write.csv(res_filter, file = paste(sub_dir, paste(savename, "csv", sep = "."), sep ="/"), quote=FALSE)
write.csv(res_subset, file = paste(sub_dir, paste(savename, "lfc1.csv", sep = "_"), sep ="/"), quote=FALSE)

postscript(file=paste(plot_dir, paste(savename, "ps", sep = "."), sep ="/"), width=900, height=600)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
png(file=paste(plot_dir, paste(savename, "png", sep = "."), sep ="/"), width=600, height=400)
plotMA(res, ylim=c(-3,3), colSig = "#c00000")
dev.off()
