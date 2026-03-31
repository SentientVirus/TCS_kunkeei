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
to_install <- c("DESeq2", "GenomeInfoDbData", "ggplot2", "stringr", "apeglm", "pheatmap", "ggrepel", "ape")
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
# 1. Load input variables from Snakemake                                      #
#=============================================================================#

flog.info("Definining input variables")
sampleFiles <- snakemake@input[["counts"]] #Input count files
sub_dir <- dirname(snakemake@output[["dif_expr"]][1]) #Directory with differential expression results
plot_dir <- dirname(snakemake@output[["pca"]][1]) #Directory with plots
pcaplot <- snakemake@output[["pca"]] #PCA plot file
distplot <- snakemake@output[["dist"]] #Distance plot file
heatplot <- snakemake@output[["heatmap"]][1] #Global heatmap file
metadat <- snakemake@input[["meta"]] #Metadata file

cw <- 18 #Heatmap cell width

# Create output directories if they don't exist
if (!file.exists(sub_dir)){
  dir.create(sub_dir)
}

if (!file.exists(plot_dir)){
  dir.create(plot_dir)
}

#=============================================================================#
# 2. Load annotations from the GenBank file                                   #
#=============================================================================#

flog.info("Loading annotations...")
# Read two GBKS, chromosome and plasmid for strain H3B1-04J
gbks <- c("OX335197", "OX335198")

# Get annotations from both GenBanks
flog.info("Read GenBank files")
my_annot <- getAnnotationsGenBank(gbks) #Retrieve gene annotations from the GenBank files
print(my_annot)

# Create a vector to store locus tags
loctags <- c()
products <- c()

# Create a new object to filter out annotations for repeat regions
new_annot <- my_annot #Create a copy of the annotations
new_annot$OX335197 <- new_annot$OX335197[new_annot$OX335197$type == "gene",] #Keep only the features annotated as genes

# Create vector with the annotations
flog.info("Store gene names in variable")
annots <- c(new_annot$OX335197$gene, new_annot$OX335198$gene) #Get only the column with annotations

# Retrieve all the locus tags and product descriptions from the GenBank file
flog.info("Store locus tags and gene descriptions in separate variables")
for (value in c(new_annot$OX335197$others, new_annot$OX335198$others)) { #Retrieve the section with locus tags and other information
  locus_tag <- str_extract(value, "(?<=locus_tag: )[A-Z0-9_]+") #Extract locus tag
  loctags <- c(loctags, locus_tag) #Append to vector
}

for (value2 in c(new_annot$OX335197$product, new_annot$OX335198$product)){ #Loop through product descriptions
  products <- c(products, value2) #Append them to the list of products
}

# Annotations to replace manually
#to_replace <- c("H3B104J_00510", "H3B104J_00520", "H3B104J_01020",
#                "H3B104J_12990", "H3B104J_13000", "H3B104J_13010",
#                "H3B104J_13020", "H3B104J_13060", "H3B104J_13070",
#                "H3B104J_13080", "H3B104J_13090", "H3B104J_13100",
#                "H3B104J_13110", "H3B104J_13120", "H3B104J_13130",
#                "H3B104J_13140", "H3B104J_13150", "H3B104J_13160",
#                "H3B104J_13170", "H3B104J_14310", "H3B104J_PKUN00040",
#                "H3B104J_PKUN00110")
#replacement <- c("adhesin_510", "GT2_520", "adhesin_1020", "GS", "BrS_13000",
#                 "BrS_13010", "GS-BrS", "wzx", "MFS_13070", "wzy", 
#                 "GT14_13090", "GT14_13100", "GT2_13110", "GT1_13120", "epsE", 
#                 "epsD", "epsC", "epsB", "epsA", "adhesin_14310", 
#                 "kukA", "nisB")

#=============================================================================#
# 3. Read counts and create metadata                                          #
#=============================================================================#

flog.info("Reading input count files")

# Assign a condition to each sample
sampleCondition <- c()
for (file in sampleFiles){
  if ((grepl("01", file, fixed = TRUE) == 1) | (grepl("02", file, fixed = TRUE) == 1)){ #If the file name contains the names of isolates 01 or 02
    sampleCondition <- c(sampleCondition, "mucoid") #Assign the mucoid condition
  }
  else{sampleCondition <- c(sampleCondition, "inhibitor")} #Else, assign the inhibitor condition
}

# Create table with metadata (sample names and condition)
sampleTable <- data.frame(sampleName = substring(sampleFiles, first = 10, last = 11),
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition) #Convert the condition variable to factor

# Read counts from files
filename <- sampleFiles[1] #Read the name of the first file
counts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid")) #Create count matrix
for (filename in sampleFiles[2:length(sampleFiles)]){ #Loop through the remaining files
  cts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid")) #Create count matrix
  counts <- cbind(counts, cts) #Append to the original count matrix
}

# Add metadata (three conditions, mucoid/inhibitor, fructose/sucrose and date)
coldata <- read.csv(metadat, sep = "\t", row.names = 1) #Read metadata file
coldata$condition2 <- factor(coldata$condition2) #Convert the three conditions to factors
coldata$condition1 <- factor(coldata$condition1)
coldata$condition3 <- factor(coldata$condition3)
coldata

# Convert the data to DESeq input format
dds_pcr <- DESeqDataSetFromMatrix(countData = counts, #Read data as a DESeq dataframe
                                  colData = coldata,
                                  design = ~ condition3 + condition1 + condition2) #Account for the three conditions, where the batch date is a fixed effect

#=============================================================================#
# 4. Run a PCA for the complete dataset                                       #
#=============================================================================#

flog.info("Running PCA")
# PCA calculations
vsd <- vst(dds_pcr, blind=FALSE) #Perform variance normalization
pcaData <- plotPCA(vsd, intgroup=c("condition1", "condition2"), returnData=TRUE) #Generate a PCA plot
percentVar <- round(100 * attr(pcaData, "percentVar")) #Get the percentage of variance of each PC

# PCA plots and formatting
p <- ggplot(pcaData, aes(PC1, PC2, color=condition2, shape=condition1)) + #Plot data
  geom_point(size=5, alpha=1, stroke=0) + #Set data point size and transparency
  scale_color_manual(values = c("S" = "#FF5733", "F" = "#0097EF")) + #Color according to carbon source
  guides(color = guide_legend(title = "Substrate", #Set plot title
                              theme = theme(legend.title = element_text(size = 24), #Set font size of legend title and text
                                            legend.text = element_text(size = 20))), 
  shape = guide_legend(title = "Morphology", #Set the data point shape by phenotype
                       theme = theme(legend.title = element_text(size = 24), #Set font size of legend title and text
                                     legend.text = element_text(size = 20)))) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + #Add x axis label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + #Add y axis label
  coord_fixed(); #Set aspect ratio
p <- p + 
  theme_minimal() + theme(axis.title = element_text(size = 30), #Set plot theme and axis title and text font sizes
                          axis.text = element_text(size = 20), 
                          panel.grid = element_blank(), #Remove panel grid
  panel.border = element_rect(fill= "transparent", linewidth = 2)); #Remove figure background

# Saving the plots to files
flog.info("Save PCA to figure")
ggsave(pcaplot[1], width = 9, height = 6);
ggsave(pcaplot[2], width = 9, height = 6);

#=============================================================================#
# 5. Create a global heatmap                                                  #
#=============================================================================#

flog.info("Creating a global heatmap")
dds_htmp <- DESeqDataSetFromMatrix(countData = counts, #DESeq input to generate a global heatmap
                                  colData = coldata,
                                  design = ~ condition1 + condition2) #Consider only isolate and carbon source

dds_htmp <- collapseReplicates(dds_htmp, #Collapse replicates
                               as.factor(str_sub(colnames(dds_htmp), 1, 5)))

global_dds <- DESeq(dds_htmp) #Run DESeq2

res_htmp <- results(global_dds) #Create results object
summary(res_htmp); #Print summary of results
res_htmp_filter <- subset(res_htmp, padj < .1); #Get only genes with an adjusted p-value < 0.1
res_htmp_filter <- subset(res_htmp_filter, baseMean >= 100); #Get only genes with an average of 100 counts or more

vsd_htmp <- vst(dds_htmp, blind=FALSE) #Perform variance normalization

htmp_select <- res_htmp_filter[order(abs(res_htmp_filter$log2FoldChange), #Sort genes by the absolute value of log2FC in decreasing order
                                    decreasing = TRUE),]
htmp_select <- rownames(htmp_select)[1:30] #Get the top 30 genes
select <- rownames(dds_htmp) %in% htmp_select #Retrieve top 30 genes

df <- as.data.frame(colData(global_dds)[,c("condition1", "condition2")]) #Create a global dataframe
colnames(df) <- c("Phenotype", "Carbon source") #Rename columns
rownames(vsd_htmp) <- str_sub(rownames(vsd_htmp), 4, -1) #Remove "AKU" from locus tags
rownames(df) <- sub("_", "", rownames(df)) #Remove underscore from row names
colnames(vsd_htmp) <- sub("_", "", colnames(vsd_htmp)) #Remove underscore from column names

get_annot <- function(res_df, loci, annotation, descriptions, prefix = ""){ #Define function inputs
  #' Function to get annotations for each gene
  #' 
  #' @description This function loops through a datafrane with a list of genes
  #' and adds annotations to each gene based on a separate file.
  #' 
  #' @param res_df output results dataframe from DESeq2. The function loops
  #' through the locus tags in the dataframe.
  #' @param loci vector. List of locus tags.
  #' @param annotation vector. List of gene annotations.
  #' @param descriptions vector. List of gene descriptions.
  #' @param prefix string. Needed to get the full locus tags, since I removed the first three characters.
  #' @usage get_annot(res_df, loci, annotations, descriptions, prefix)
  #' @return double with gene names and gene descriptions as columns.
  
  gene_names <- c() #Initialize vectors
  prod <- c()
  print(vsd_htmp)
  for (rown in rownames(res_df)){ #Loop through rows in input dataframe
    add_names <- annotation[paste(prefix, rown, sep = "") == loci][1] #Retrieve the annotation for the gene
    prod <- c(prod, descriptions[rown == loci][1]) #Append the product description of the gene
    if (!is.na(add_names)){ #If the add_names vector is not NA (not a hypothetical protein)
      gene_names <- c(gene_names, add_names) #Append it to the vector with gene names
    } else {
      gene_names <- c(gene_names, "-") #Otherwise, append - to the vector with gene names
    }
  }
  return(cbind(gene_names, prod)) #Return gene names and product descriptions
}

# Get annotations for each locus tag
vsd_info <- get_annot(vsd_htmp, loctags, annots, products, "AKU")

flog.info("Add annotations to heatmap")
# Replace the locus tags to the manual annotations
#rownames(vsd_htmp)[which(rownames(vsd_htmp) %in% to_replace)] <- replacement

#Change the row names of the genes that are not hypothetical proteins to their annotations
rownames(vsd_htmp) <- ifelse(!vsd_info[, 1] == "-", vsd_info[, 1], rownames(vsd_htmp))

flog.info("Save global heatmap")
png(file = heatplot, width = 600, height = 400); #Create PNG file
pheatmap(assay(vsd_htmp)[select,], cluster_rows=FALSE, show_rownames=TRUE, #Plot heatmap to file
         cluster_cols=TRUE, annotation_col=df, cellwidth=40)
invisible(dev.off()) #Restart canvas

postscript(file = snakemake@output[["heatmap"]][2], width = 600, height = 400); #Create postcript file
pheatmap(assay(vsd_htmp)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cellwidth=40)
invisible(dev.off())

sampleDists <- dist(t(assay(vsd_htmp)))

#=============================================================================#
# 6. Create a distance plot                                                   #
#=============================================================================#

flog.info("Creating a distance plot")
sampleDistMatrix <- as.matrix(sampleDists) #Create a matrix of distances
colors <- colorRampPalette( rev(brewer.pal(9, "BuPu")) )(255) #Assign a color palette

png(file = distplot[1], width = 600, height = 400); #Create a PNG file
postscript(file = distplot[2], width = 600, height = 400); #And a postcript file
pheatmap(sampleDistMatrix, #Plot heatmap of distances
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
invisible(dev.off()) #Restart canvas

#=============================================================================#
# 7. Differential expression analyses                                         #
#=============================================================================#

# Define lists with information to loop over
comparisons <- c("Smucoid_vs_Sinhibitor", "Fmucoid_vs_Finhibitor",
                 "Smucoid_vs_Fmucoid", "Sinhibitor_vs_Finhibitor")

coefs <- c("condition1_mucoid_vs_inhibitor", "condition1_mucoid_vs_inhibitor",
           "condition2_S_vs_F", "condition2_S_vs_F")

my_condition <- c("condition1", "condition1", "condition2", "condition2")
labels <- c("Phenotype", "Phenotype", "Carbon source", "Carbon source")
  
grep_patterns <- c("_S_", "_F_", "I01|I02", "I09|I10")

# Loop to generate results for each comparison
for (i in 1:length(comparisons)){ #Loop through comparisons
  savename <- comparisons[i] #Get the base name of the output file
  flog.info(paste("Running DESeq2 for comparison", savename))
  
  # Get data for the differential expression analysis
  if (grepl("mucoid", savename) & grepl("inhibitor", savename)){ #If it is a mucoid vs inhibitor comparison
    if (grepl("S", savename)){ #And the carbon source is sucrose
      flog.info("Muc vs Inh, S+")
      coldata2 <- coldata[coldata$condition2 == "S",] #Set the carbon source condition to keep
    }
    else { #Same, but without sucrose
      flog.info("Muc vs Inh, S-")
      coldata2 <- coldata[coldata$condition2 == "F",]
    }
  # Retrieve data for the comparison of interest only
  coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
  counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
  coldata2$condition3 = droplevels(coldata2$condition3) #The batch date is disregarded
    
  # Create DESeq2 input
  dds <- DESeqDataSetFromMatrix(countData = counts2, #Create input dataframe for DESeq2
                                  colData = coldata2,
                                  design = ~ condition3 + condition1)
  }
  else { #If the comparison is based on carbon source
    if(grepl("mucoid", savename)){ #If the phenotype is mucoid
      flog.info("S+ vs S-, Muc") #Set title
      coldata2 <- coldata[coldata$condition1 == "mucoid",] #Set condition to retrieve
    }
    else { #Same, but for the other phenotype
      flog.info("S+ vs S-, Inh")
      coldata2 <- coldata[coldata$condition1 == "inhibitor",]
    }
  # Retrieve data for the condition of interest
  coldata2 <- coldata2[rownames(coldata2) %in% colnames(counts),]
  counts2 <- counts[, colnames(counts) %in% rownames(coldata2)]
  coldata2$condition3 = droplevels(coldata2$condition3)
    
  # Create DESeq2 input 
  dds <- DESeqDataSetFromMatrix(countData = counts2, colData = coldata2,
                                design = ~ condition3 + condition2)
}

# Run the differential expression analysis
flog.info("Run DESeq2")
dds <- DESeq(dds) #Run DESeq2

# Get results with and without lfc threshold filters
flog.info("Add lfc shrinkage and threshold")
res <- lfcShrink(dds, coef=coefs[i], type="apeglm") #Get results applying log2FC shrinkage
res_subset <- results(dds, lfcThreshold=1); #Filter out results with log2FC < 1
summary(res); #Print a summary of the results

flog.info("Filter by p-value")
res_subset <- subset(res_subset, padj < .1); #Filter out results with an adjusted p-value > 0.1
res_filter <- subset(res, padj < .1);
res_comp_htmp <- subset(res_filter, baseMean >= 100); #Filter out results with mean counts < 100

# Create a heatmap for this comparison
flog.info(paste("Create a heatmap for comparison", savename))
res_select <- res_comp_htmp[order(abs(res_comp_htmp$log2FoldChange), decreasing = TRUE),] #Sort results by the absolute log2FC value
res_select <- rownames(res_select)[1:20] #Get the 20 genes with the highest log2FC in absolute value
select <- rownames(dds) %in% res_select #Filter out the remaining genes

df <- as.data.frame(colData(dds)[,c(my_condition[i], "condition3")]) #Create a dataframe from the DESeq output where the conditions are the column names
colnames(df) <- c(labels[i], "Batch") #Change column names
vsd_subset <- vsd[, sapply(colnames(vsd), grepl, pattern = grep_patterns[i])] #Perform variance transformation on the data of interest
rownames(vsd_subset) <- str_sub(rownames(vsd_subset), 4, -1) #Remove "AKU" from locus tags

# Get annotations for each locus tag
vsd_info <- get_annot(vsd_subset, loctags, annots, products, "AKU") #Retrieve annotations and product descriptions
# Replace locus tag for manual annotations
#rownames(vsd_subset)[which(rownames(vsd_subset) %in% to_replace)] <- replacement
rownames(vsd_subset) <- ifelse(!vsd_info[, 1] == "-", vsd_info[, 1], rownames(vsd_subset)) #Change the label of the genes from the locus tag to the annotations when the protein is not hypothetical

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
write.csv(res_filter, file = snakemake@output[["dif_expr"]][2*i-1], quote=FALSE); #Write DESeq2 statistics to CSV
write.csv(res_subset, file = snakemake@output[["dif_expr"]][2*i], quote=FALSE);

# Section to annotate output
flog.info(paste("Add annotations to", savename))
to_annotate <- list(res_filter, res_subset) #Set output to annotate

# Implement function to add annotations
k <- 2*i-1 #Index to get the right output file
for (my_res in to_annotate){ #Loop through outputs to annotate
  res_info <- get_annot(my_res, loctags, annots, products) #Retrieve annotations for each output
  my_res$gene_names <- res_info[, 1] #Add new column with annotations
  my_res$products <- res_info[ , 2] #Add new column with product descriptions
  write.table(my_res, file = snakemake@output[["DE_annot"]][k], quote=FALSE, sep = "\t", col.names = NA); #Write annotated results to file
  k <- k + 1
  }

# Generate MA plots (png and postscript)
png(file = snakemake@output[["plots"]][2*i-1], width = 600, height = 400);
plotMA(res, ylim = c(-3,3), colSig = "#c00000"); #Plot a plotMA where significant dots are displayed in red
invisible(dev.off())
postscript(file = snakemake@output[["plots"]][2*i], width = 900, height = 600);
plotMA(res, ylim = c(-3,3), colSig = "#c00000");
invisible(dev.off())}