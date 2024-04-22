#=============================================================================#
# 0. Install packages if needed                                               #
#=============================================================================#
list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org');

to_install <- c("ggplot2", "ggrepel", "EnhancedVolcano")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package);
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#
#library("DESeq2")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
#library("pheatmap")
#library("RColorBrewer")
#library("stringr")
#library(futile.logger)

#=============================================================================#
# 1. Load dataframe with DESeq2 output + annotations                          #
#=============================================================================#

# Script to generate a Volcano plot using ggplot, following this tutorial:
# https://biostatsquid.com/volcano-plots-r-tutorial/

#cols <- c("character",  rep("numeric", 10), rep("numeric", 10), rep("numeric", 10), rep("numeric", 10), rep("numeric", 10), "character", "character")
res <- read.table("../results/DE/Smucoid_vs_Sinhibitor_annotated.tsv", sep = "\t", 
                 numerals = "no.loss", header = TRUE, row.names = 1, quote = "",
                 stringsAsFactors = FALSE)


#=============================================================================#
# 2. Volcano plot with ggplot and standard plotting                           #
#=============================================================================#
yax <- -log10(as.numeric(res$padj))
xax <- as.numeric(res$log2FoldChange)
yval <- yax
yval[!is.finite(yval)] <- 300

# Variable that I use to change the shape of the points
sign_shape <- sign(xax)
sign_shape[sign_shape == 1] <- "Mucoid"
sign_shape[sign_shape == -1] <- "Inhibitor"
sign_shape <- as.factor(sign_shape)
#sign_color <- sign(xax)
#sign_color[sign_color > 0] <- 0
#sign_color[abs(xax) < 0.5] <- 1
#sign_color[yval < -log10(0.05)] <- 2
#sign_color <- as.factor(sign_color)

# Variable that I use to annotate genes in the plots
vollabels <- res$gene_name # Get gene names
# Label genes of interest manually
vollabels["AKUH3B104J_00510" == rownames(res)] <- "adhesin_510"  #expression("adhesin"["510"])
vollabels["AKUH3B104J_00520" == rownames(res)] <- "adhesin_520"  #expression("adhesin"["PNAG"])
vollabels["AKUH3B104J_01020" == rownames(res)] <- "adhesin_1020"  #expression("adhesin"["1020"])
vollabels["AKUH3B104J_14310" == rownames(res)] <-  "adhesin_14310" #expression("adhesin"["14310"])
vollabels["AKUH3B104J_PKUN00040" == rownames(res)] <- "kukA"
# Add labels only in genes that are not hypothetical and have x > 3 and y > 80
vollabels[(abs(xax) < 3) & (yax < 80)] <- ""
vollabels[vollabels == "-"] <- ""
# Make labels bold and italic
italic_labels <- vollabels
italic_labels[!italic_labels == ""] <- paste0("bolditalic('", vollabels[!vollabels == ""],"')")
select_labs <- italic_labels[!italic_labels == ""]

# Variable that I use to convert some variables in res to numeric
new_res <- res
new_res$log2FoldChange <- as.numeric(new_res$log2FoldChange)
new_res$padj <- as.numeric(new_res$padj)
new_res$pvalue <- as.numeric(new_res$pvalue)

# Variable that I use to store the colors
keyvals.col <- c()
# Coloring depending on x (different colors if it's < 0.5, > 3 or between) and y
# (different colors depending on if the padj is <1e-5 or >1e-5)
keyvals.col <- ifelse(
  abs(new_res$log2FoldChange) < 0.5 & new_res$padj < 1e-5, "darkblue",
  ifelse(abs(new_res$log2FoldChange) > 0.5 & abs(new_res$log2FoldChange) < 3 & new_res$padj < 1e-5, "cyan4",
         ifelse(abs(new_res$log2FoldChange) > 3 & new_res$padj < 1e-5, "aquamarine3",
                ifelse(abs(new_res$log2FoldChange) < 0.5 & new_res$padj > 1e-5, "bisque4",
                       ifelse(abs(new_res$log2FoldChange) > 3 & new_res$padj > 1e-5, "bisque2",
                              "bisque3")))))

# Labels for each color
keyvals.col[is.na(keyvals.col)] <- "red"
names(keyvals.col)[keyvals.col == "aquamarine3"] <- expression(italic("p"["adj"])*" < 10"^-5*", Log"[2]*italic("FC")*" > 3")
names(keyvals.col)[keyvals.col == "cyan4"] <- expression(italic("p"["adj"])*" < 10"^-5*", 0.5 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "darkblue"] <- expression(italic("p"["adj"])*" < 10"^-5*", Log"[2]*italic("FC")*" < 0.5")
names(keyvals.col)[keyvals.col == "bisque2"] <- expression(italic("p"["adj"])*" > 10"^-5*", Log"[2]*italic("FC")*" < 0.5")
names(keyvals.col)[keyvals.col == "bisque3"] <- expression(italic("p"["adj"])*" > 10"^-5*", 0.5 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "bisque4"] <- expression(italic("p"["adj"])*" > 10"^-5*", Log"[2]*italic("FC")*" < 0.5")

color_values <- unique(keyvals.col)
color_values <- color_values[order(color_values)] #c("aquamarine3", "bisque3", "bisque4", "cyan4", "darkblue")
labels <- c(unique(names(keyvals.col)[keyvals.col == color_values[1]]), 
            unique(names(keyvals.col)[keyvals.col == color_values[2]]), 
            unique(names(keyvals.col)[keyvals.col == color_values[3]]), 
            unique(names(keyvals.col)[keyvals.col == color_values[4]]), 
            unique(names(keyvals.col)[keyvals.col == color_values[5]]))

size_vector <- abs(res$log2FoldChange)

volcanoplot <- ggplot(data = res, aes(x = log2FoldChange, y = yval, col = keyvals.col, label = vollabels)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = "dashed") + # Add dashed line to show log2FC < 0.5
  geom_hline(yintercept = -log10(1e-5), col = "gray", linetype = "dashed") + # Add dashed line for p-value > 0.05
  geom_point(aes(size = size_vector, shape = sign_shape, fill = keyvals.col), alpha = 0.6, stroke = 0.5, color = "brown1") + 
  scale_shape_manual(values = c(21, 24)) + #, labels = c("Inhibitor", "Mucoid")) +
  scale_size(range = c(0.5, 3)) + #, guide = "none") + 
  theme_light() + # Set point size and overall graph appearance
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # Remove grid
  scale_fill_manual(values = color_values, # Set the colors of up/downregulated points
      labels = parse(text = labels), aesthetics = c("colour", "fill")) + # Set the color labels
  #scale_fill_manual(values = color_values) +
  guides(fill = guide_legend(override.aes = list(size = 3)), # Increase the size of legend points
         shape = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(ylim = c(0, 320), xlim = c(-8, 8)) + # Set plot limits
  labs(fill = "Differential expression", shape = "Morphology", size = expression("Log"[2]*italic("Fold Change")), # Legend title
      x = expression("Log"[2]*italic("Fold Change")), y = expression("-Log"[10]*italic("p"["adj"]))) +
  scale_x_continuous(breaks = seq(-8, 8, 2)) + # Customize ticks in the x axis
  scale_y_continuous(breaks = seq(0, 320, 40)) + # Customize ticks in y axis
  ggtitle("Mucoid vs Inhibitor", subtitle = "+ sucrose") + # Plot title
  geom_label_repel(max.overlaps = Inf, show_guide = FALSE, color = "black", 
      size = 3, box.padding = 0.4, fontface = "bold.italic") #, label.size = NA, fill = NA) # To show all labels

volcanoplot <- volcanoplot + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
volcanoplot
# Saving the plots to files
ggsave("../try.png", width = 9, height = 6)


#=============================================================================#
# 3. Volcano plot with specific package for this                              #
#=============================================================================#
new_res <- res
new_res$log2FoldChange <- as.numeric(new_res$log2FoldChange)
new_res$padj <- as.numeric(new_res$padj)
new_res$pvalue <- as.numeric(new_res$pvalue)
italic_labels <- vollabels
italic_labels[!italic_labels == ""] <- paste0("bolditalic('", vollabels[!vollabels == ""],"')")
select_labs <- italic_labels[!italic_labels == ""]

keyvals.colour <- ifelse(
  abs(new_res$log2FoldChange) < 0.5 & new_res$padj < 1e-5, "darkblue",
  ifelse(abs(new_res$log2FoldChange) > 0.5 & abs(new_res$log2FoldChange) < 3 & new_res$padj < 1e-5, "cyan4",
  ifelse(abs(new_res$log2FoldChange) > 3 & new_res$padj < 1e-5, "aquamarine3",
  ifelse(abs(new_res$log2FoldChange) < 0.5 & new_res$padj > 1e-5, "bisque4",
  ifelse(abs(new_res$log2FoldChange) > 3 & new_res$padj > 1e-5, "bisque2",
         "bisque3")))))
keyvals.colour[is.na(keyvals.colour)] <- "red"
names(keyvals.colour)[keyvals.colour == "aquamarine3"] <- "Significant, high" #expression("p"["adj"]) #*" < 1e"^(-5)*", Log"[2]~"FC > 3")
names(keyvals.colour)[keyvals.colour == "cyan4"] <- "Significant, mid"
names(keyvals.colour)[keyvals.colour == "darkblue"] <- "Significant, low"
names(keyvals.colour)[keyvals.colour == "bisque2"] <- "Not significant, high"
names(keyvals.colour)[keyvals.colour == "bisque3"] <- "Not significant, mid"
names(keyvals.colour)[keyvals.colour == "bisque4"] <- "Not significant, low"

enhanced_plot <- EnhancedVolcano(new_res, lab = italic_labels, x = "log2FoldChange",
    y = "padj", pCutoff = 1e-5, FCcutoff = 0.5,  
    hline = 10e-6, hlineCol = "grey", hlineType = "longdash",
    vline = c(-0.5, 0.5), vlineCol = "grey", vlineType = "longdash",
    legendPosition = 'right', legendLabSize = 12, legendIconSize = 4.0,
    drawConnectors = TRUE, widthConnectors = 1, colConnectors = 'black',
    labSize = 3.0, labCol = 'black', labFace = 'bold', boxedLabels = TRUE,
    parseLabels = TRUE, selectLab = select_labs, max.overlaps = Inf,
    border = "full", gridlines.minor = FALSE, gridlines.major = FALSE,
    title = "Mucoid vs inhibitor", subtitle = "+ Sucrose", 
    colAlpha = 0.7, colCustom = keyvals.colour, legendLabels = parse(names(keyvals.colour)), 
    pointSize = c(ifelse(abs(new_res$log2FoldChange) > 2.5 | new_res$padj <= 1e-80, 2, 1)))
enhanced_plot <- enhanced_plot +
    ggplot2::coord_cartesian(xlim=c(-8, 8)) +
    ggplot2::scale_x_continuous(breaks=seq(-8, 8, 2))
enhanced_plot
ggsave("../try2.png", width = 9, height = 6)