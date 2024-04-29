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

to_install <- c("DESeq2", "ggplot2", "stringr", "apeglm", "pheatmap", "ggrepel", "ape")
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
library("ape")
library(futile.logger)

#=============================================================================#
# 0. Logging                                                                  #
#=============================================================================#
# Create a logger that will be saved to a file
flog.logger("saturation", TRACE, appender=appender.file(snakemake@log[[1]]))

flog.info("R script to run a differential expression analysis")

#=============================================================================#
# Provisional section to load GenBank file                                    #
#=============================================================================#
# Read two GBKS, chromosome and plasmid for strain H3B1-04J
gbks <- c("OX335197", "OX335198")

# Get annotations from both GenBanks
my_annot <- getAnnotationsGenBank(gbks)

# Create a vector to store locus tags
loctags <- c()
products <- c()

# Create a new object to filter out annotations for repeat regions
new_annot <- my_annot
new_annot$OX335197 <- new_annot$OX335197[new_annot$OX335197$type == "gene",]

# Create vector with the annotations
annots <- c(new_annot$OX335197$gene, new_annot$OX335198$gene)

# Retrieve all the locus tags from the GenBank information
for (value in c(new_annot$OX335197$others, new_annot$OX335198$others)) {
  locus_tag <- str_extract(value, "(?<=locus_tag: )[A-Z0-9_]+")
  loctags <- c(loctags, locus_tag)
}

for (value2 in c(new_annot$OX335197$product, new_annot$OX335198$product)){
  products <- c(products, value2)
}

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
  geom_point(size=8, alpha=1, stroke=0) + scale_color_manual(values = c("S" = "#FF5733", "F" = "#0097EF")) +
  guides(color = guide_legend(title = "Substrate", theme = theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20))), 
  shape = guide_legend(title = "Morphology", theme = theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20)))) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed();
p <- p + 
  #stat_ellipse(geom="polygon", aes(fill = pcaData$condition2), 
  #alpha = 0.2, show.legend = FALSE, level = 0.95) + 
  scale_fill_manual(values = c("S" = "#FF5733", "F" = "#0097EF")) +
  theme_minimal() + theme(axis.title = element_text(size = 30), axis.text = element_text(size = 20), panel.grid = element_blank(), 
  panel.border = element_rect(fill= "transparent", size = 2));

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

# Get annotations for each locus tag
my_rownames <- c()
for (rown in rownames(vsd_htmp)){
  comp <- annots[paste("AKU", rown, sep = "") == loctags][1]
  if (!is.na(comp)){
    my_rownames <- c(my_rownames, comp)
  } else {
    my_rownames <- c(my_rownames, rown)
    }
}

flog.info("Add annotations to heatmap")
rownames(vsd_htmp) <- my_rownames

# Replace the names of certain rows
to_replace <- c("H3B104J_00510", "H3B104J_00520", "H3B104J_01020", "H3B104J_14310", "H3B104J_PKUN00040")
replacement <- c("adhesin_510", "adhesin_520", "adhesin_1020", "adhesin_14310", "kukA")

# Implement replacement
rownames(vsd_htmp)[which(rownames(vsd_htmp) %in% to_replace)] <- replacement

flog.info("Save global heatmap")
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
comparisons <- c("Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor",
                 "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor")

coefs <- c("condition1_mucoid_vs_inhibitor", "condition1_mucoid_vs_inhibitor",
           "condition2_S_vs_F", "condition2_S_vs_F")

my_condition <- c("condition1", "condition1", "condition2", "condition2")
labels <- c("Phenotype", "Phenotype", "Carbon source", "Carbon source")
  
grep_patterns <- c("_S_", "_F_", "I01|I02", "I09|I10")

for (i in 1:length(comparisons)){
  # Run a differential expression analysis for each comparison
  savename <- comparisons[i]
  flog.info(paste("Running DESeq2 for comparison", savename))
  
  # Get data for the differential expression analysis
  if (grepl("mucoid", savename) & grepl("inhibitor", savename)){
    if (grepl("S", savename)){
      flog.info("Muc vs Inh, S+")
      coldata2 <- coldata[coldata$condition2 == "S",]
    }
    else {
      flog.info("Muc vs Inh, S-")
      coldata2 <- coldata[coldata$condition2 == "F",]
    }
  coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
  counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
  coldata2$condition3 = droplevels(coldata2$condition3)
    
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts2,
                                  colData = coldata2,
                                  design = ~ condition3 + condition1)
  }
  else {
    if(grepl("mucoid", savename)){
      flog.info("S+ vs S-, Muc")
      coldata2 <- coldata[coldata$condition1 == "mucoid",]
    }
    else {
      flog.info("S+ vs S-, Inh")
      coldata2 <- coldata[coldata$condition1 == "inhibitor",]
    }
  coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
  counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
  coldata2$condition3 = droplevels(coldata2$condition3)
    
  dds <- DESeqDataSetFromMatrix(countData = counts2, colData = coldata2,
                                design = ~ condition3 + condition2)
}

# Run the differential expression analysis
flog.info("Run DESeq2")
dds <- DESeq(dds)

# Get results with and without lfc threshold filters
flog.info("Add lfc shrinkage and threshold")
res <- lfcShrink(dds, coef=coefs[i], type="apeglm")
res_subset <- results(dds, lfcThreshold=1);
summary(res);

flog.info("Filter by p-value")
res_subset <- subset(res_subset, padj < .1);
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100);

# Create a heatmap for this comparison
flog.info(paste("Create a heatmap for comparison", savename))
res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),]
res_select <- rownames(res_select)[1:20]
select <- rownames(dds) %in% res_select

df <- as.data.frame(colData(dds)[,c(my_condition[i], "condition3")])
colnames(df) <- c(labels[i], "Batch")
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = grep_patterns[i])]
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1)

# Get annotations for each locus tag
my_rownames <- c()
for (rown in rownames(vsd_subset)){
  comp <- annots[paste("AKU", rown, sep = "") == loctags][1]
  if (!is.na(comp)){
    my_rownames <- c(my_rownames, comp)
  } else {
    my_rownames <- c(my_rownames, rown)
  }
}
rownames(vsd_subset) <- my_rownames
to_replace <- c("H3B104J_00510", "H3B104J_00520", "H3B104J_01020", "H3B104J_13020", "H3B104J_14310", "H3B104J_PKUN00040")
replacement <- c("adhesin_510", "adhesin_520", "adhesin_1020", "GS2_BRS", "adhesin_14310", "kukA")
rownames(vsd_subset)[which(rownames(vsd_subset) %in% to_replace)] <- replacement

# Save heatmap to files
png(file = snakemake@output[["heatmap"]][2*i+1], width = 600, height = 400) 
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

postscript(file = snakemake@output[["heatmap"]][2*(i+1)], width = 600, height = 400)
pheatmap(assay(vsd_subset)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=cw)
invisible(dev.off())

# Save unnanotated output to file
flog.info(paste("Save DE output (without annotations) to file for comparison", savename))
write.csv(res_filter, file = snakemake@output[["dif_expr"]][2*i-1], quote=FALSE);
write.csv(res_subset, file = snakemake@output[["dif_expr"]][2*i], quote=FALSE);

# Add annotations to output
flog.info(paste("Add annotations to", savename))
to_annotate <- list(res_filter, res_subset)

get_annot <- function(res_df, loci, annotation, descriptions){
  gene_names <- c()
  prod <- c()
  for (rown in rownames(res_df)){
    add_names <- annotation[rown == loci][1]
    prod <- c(prod, descriptions[rown == loci][1])
    if (!is.na(add_names)){
      gene_names <- c(gene_names, add_names)
    } else {
      gene_names <- c(gene_names, "-")
    }
  }
  return(cbind(gene_names, prod))
}

resf_info <- get_annot(res_filter, loctags, annots, products)
res_filter$gene_names <- resf_info[, 1]
res_filter$products <- resf_info[, 2]
ress_info <- get_annot(res_subset, loctags, annots, products)
res_subset$gene_names <- ress_info[, 1]
res_subset$products <- ress_info[, 2]

#new_cols <- c("locus_tag", colnames(res_filter))
#print(new_cols)

# Save annotated results to tsv files
write.table(res_filter, file = snakemake@output[["DE_annot"]][2*i-1], quote=FALSE, sep = "\t", col.names = NA);
write.table(res_subset, file = snakemake@output[["DE_annot"]][2*i], quote=FALSE, sep = "\t", col.names = NA);

# Generate MA plots (png and postscript)
png(file = snakemake@output[["plots"]][2*i-1], width = 600, height = 400); #OBS!
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())
postscript(file = snakemake@output[["plots"]][2*i], width = 900, height = 600); #OBS!
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())}