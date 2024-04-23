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
# Define inputs
input_files <- c("../results/DE/Smucoid_vs_Sinhibitor_annotated.tsv",
                 "../results/DE/Fmucoid_vs_Finhibitor_annotated.tsv",
                 "../results/DE/Smucoid_vs_Fmucoid_annotated.tsv",
                 "../results/DE/Sinhibitor_vs_Finhibitor_annotated.tsv")

# Define vectors with formatting
titles <- rep(c("Mucoid vs Inhibitor", " + Sucrose vs - Sucrose"), each = 2)
subtitles <- c("+ sucrose", "- sucrose", "Mucoid", "Inhibitor")
all_shapes <- rep(data.frame(c(21, 23), c(24, 25)), each = 2)
comparison <- rep(data.frame(c("Inh", "Muc"), c("+", "-")), each = 2)
size_breaks_all <- rep(data.frame(c(2, 4, 6), c(1, 2, 3)), each = 2)
scales <- rep(data.frame(c(0.5, 3), c(0.25, 1.5)), each = 2)
shape_names <- rep(c("Morphology", "Sucrose"), each = 2)

# Loop through i to get res and plot
for (i in 1:4){

# Read the annotated DESeq2 output
res <- read.table(input_files[i], sep = "\t", 
                 numerals = "no.loss", header = TRUE, row.names = 1, quote = "",
                 stringsAsFactors = FALSE)

# Convert numeric columns to numeric (character by default)
j <- 1:5
res[, j] <- apply(res[, j], 2, function(x) as.numeric(x))

# Get shapes and titles for the plot
shapes <- as.vector(unlist(all_shapes[i]))
size_breaks <- as.vector(unlist(size_breaks_all[i]))
scale <- as.vector(unlist(scales[i]))
size_labels <- apply(expand.grid(size_breaks, as.vector(unlist(comparison[i]))), 1, paste, collapse=", ") # Set different values for different comparisons
title <- titles[i]
subtitle <- subtitles[i]
shape_label <- shape_names[i]

#=============================================================================#
# 2. Volcano plot with ggplot and standard plotting                           #
#=============================================================================#
# Set the variable 
yax <- -log10(res$padj)
xax <- res$log2FoldChange
yval <- yax
ymax <- max(yval[is.finite(yval)])
yval[!is.finite(yval)] <- ymax
xmax <- max(abs(xax))


# Variable used to change the shape of the points
sign_shape <- sign(xax)
# TO DO: Adapt these lines so they are all set to +-sucrose (and to the + and - signs) for the other comparisons
sign_shape[sign_shape == 1] <- shapes[2]
sign_shape[sign_shape == -1] <- shapes[1]
sign_shape <- as.factor(sign_shape)

# Variable used to annotate genes in the plots
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

# Variable that I use to store the colors
keyvals.col <- c()
# Coloring depending on x (different colors if it's < 0.5, > 3 or between) and y
# (different colors depending on if the padj is <1e-5 or >1e-5)
keyvals.col <- ifelse(
  abs(res$log2FoldChange) < 0.5 & res$padj < 1e-5, "#00008B",
  ifelse(abs(res$log2FoldChange) > 0.5 & abs(res$log2FoldChange) < 3 & res$padj < 1e-5, "#008B8B",
         ifelse(abs(res$log2FoldChange) > 3 & res$padj < 1e-5, "#66CDAA",
                ifelse(abs(res$log2FoldChange) < 0.5 & res$padj > 1e-5, "#8B7D6B",
                       ifelse(abs(res$log2FoldChange) > 3 & res$padj > 1e-5, "#EED5B7",
                              "#CDB79E")))))

# Labels for each color
keyvals.col[is.na(keyvals.col)] <- "red"
names(keyvals.col)[keyvals.col == "#66CDAA"] <- expression(italic("p"["adj"])*" < 10"^-5*", Log"[2]*italic("FC")*" > 3")
names(keyvals.col)[keyvals.col == "#008B8B"] <- expression(italic("p"["adj"])*" < 10"^-5*", 0.5 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "#00008B"] <- expression(italic("p"["adj"])*" < 10"^-5*", Log"[2]*italic("FC")*" < 0.5")
names(keyvals.col)[keyvals.col == "#EED5B7"] <- expression(italic("p"["adj"])*" > 10"^-5*", Log"[2]*italic("FC")*" < 0.5")
names(keyvals.col)[keyvals.col == "#CDB79E"] <- expression(italic("p"["adj"])*" > 10"^-5*", 0.5 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "#8B7D6B"] <- expression(italic("p"["adj"])*" > 10"^-5*", Log"[2]*italic("FC")*" < 0.5")

# Assign colors to data points
color_values <- unique(keyvals.col)
color_values <- color_values[order(color_values)] # Colors must be ordered alphabetically to match the plot
labels <- unlist(lapply(color_values, function(color) {
  unique(names(keyvals.col)[keyvals.col == color])}))
  
#labels <- c(unique(names(keyvals.col)[keyvals.col == color_values[1]]), # Add the right label to each color
#            unique(names(keyvals.col)[keyvals.col == color_values[2]]), 
#            unique(names(keyvals.col)[keyvals.col == color_values[3]]), 
#            unique(names(keyvals.col)[keyvals.col == color_values[4]]), 
#            unique(names(keyvals.col)[keyvals.col == color_values[5]]))

# Vector to scale the size of data points
size_vector <- abs(res$log2FoldChange)

# Code to generate the Volcano plot
volcanoplot <- ggplot(data = res, aes(x = log2FoldChange, y = yval, col = keyvals.col, label = vollabels)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = "dashed") + # Add dashed line to show log2FC < 0.5
  geom_hline(yintercept = -log10(1e-5), col = "gray", linetype = "dashed") + # Add dashed line for p-value > 1e-5
  geom_point(aes(size = size_vector, shape = sign_shape, fill = keyvals.col), alpha = 0.6, stroke = 0.5, color = "darkorchid") + # Line the size, shape, color (fill + border) and stroke of the points
  geom_label_repel(max.overlaps = Inf, show_guide = FALSE, color = "black", 
        size = 3, box.padding = 0.4, fontface = "bold.italic") + # Add label boxes
  scale_shape_manual(values = shapes, guide = "none") + # TO DO: Make a shape vector to use different shapes depending on the comparison
  #scale_size(range = c(0.5, 3)) + # Scale point size so that it is not too big
  theme_light() + # Set point size and overall graph appearance
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # Remove grid
  scale_size_continuous(breaks = rep(size_breaks, times = 2), labels = size_labels, range = scale) +
  scale_fill_manual(values = color_values, # Set the colors of up/downregulated points
      labels = parse(text = labels)) + # Set the color labels
  guides(fill = guide_legend(override.aes = list(order = 3, size = 3, shape = 21, color = "darkorchid")), # Increase the size of legend points, change shape to add border color
       #shape = guide_legend(override.aes = list(order = 2, size = 3, color = "darkorchid")), # Same as above, but for the shape legend
       size = guide_legend(nrow = 3, ncol = 2, bycol = TRUE, 
       override.aes = list(shape = rep(shapes, each = 3), order = 1, color = "darkorchid"))) + # Change the color of the size legend to gray, TO DO: Try to add triangles as well to the legend
  labs(fill = "Differential expression", shape = shape_label, 
       size = expression("Log"[2]*italic("Fold Change")), # Legend title
       x = expression("Log"[2]*italic("Fold Change")), # Title of main axes
       y = expression("-Log"[10]*italic("p"["adj"])), color = FALSE) +
  coord_cartesian(ylim = c(0, 10*round(ymax/10)), xlim = c(-(round(xmax) + 1), round(xmax) + 1)) + # Set axis limits
  scale_x_continuous(breaks = seq(-(round(xmax) + round(xmax) + 1), 8, 1)) + # Customize ticks in the x axis
  scale_y_continuous(breaks = seq(0, 10*round(ymax/10), 25)) + # Customize ticks in y axis
  ggtitle(title, subtitle = subtitle) + # Plot title
  theme(plot.title = element_text(hjust = 0.5), # Center title
       plot.subtitle = element_text(hjust = 0.5)) # Center subtitle

# Show plot in console
volcanoplot

# Saving the plots to files, TO DO: save each plot a different file, maybe using a function
ggsave(paste("../", paste(i, ".png", sep = ""), sep = ""), width = 9, height = 6)}

