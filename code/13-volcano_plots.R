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

to_install <- c("ggplot2", "ggrepel")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package);
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#

library("ggplot2")
library("ggrepel")
library(futile.logger)

#=============================================================================#
# 0. Logging                                                                  #
#=============================================================================#

# Create a logger that will be saved to a file
flog.logger("saturation", TRACE, appender=appender.file(snakemake@log[[1]]))

flog.info("R script to generate volcano plots")

#=============================================================================#
# 1. Load dataframe with DESeq2 output + annotations                          #
#=============================================================================#

# Define inputs
flog.info("Definining input variables")
input_files <- snakemake@input[["annotated_expr"]] #Annotated DESeq2 results
output_files <- snakemake@output[["volcano"]] #Volcano plot files

# Define vectors with formatting
titles <- rep(c("Aggregating vs Mucoid", " - Sucrose vs + Sucrose"), each = 2) # Plot titles
subtitles <- c("+ sucrose", "- sucrose", "Mucoid", "Aggregating") # Plot subtitles
all_shapes <- rep(data.frame(c(21, 23), c(25, 24)), each = 2) # Desired point shape for the plot
comparison <- rep(data.frame(c("Agg", "Muc"), c("-", "+")), each = 2) # Comparisons to be plotted
scale_val <- c(0, 3) # Variable to scale points

pval_filter <- 0.1 #p-adjusted threshold
lfc_filter <- 1 #Log2FC threshold

# Loop through i to get res and plot
for (i in 1:4){
  
# Read the annotated DESeq2 output
flog.info(paste("Reading file ", input_files[i])) 
res <- read.table(input_files[i], sep = "\t", 
                 numerals = "no.loss", header = TRUE, row.names = 1, quote = "",
                 stringsAsFactors = FALSE)

# Convert numeric columns to numeric (character by default)
j <- 1:5
res[, j] <- apply(res[, j], 2, function(x) as.numeric(x))

# Get shapes and titles for the plot
shapes <- as.vector(unlist(all_shapes[i]))
title <- titles[i]
subtitle <- subtitles[i]

#=============================================================================#
# 2. Volcano plot with ggplot and standard plotting                           #
#=============================================================================#

# Set the x and y variables 
flog.info("Calculate x and y") 
yax <- -log10(res$padj)
xax <- res$log2FoldChange

# Get maximum and set infinite values (padj = 0) to maximum
yval <- yax
ymax <- max(yval[is.finite(yval)])
yval[!is.finite(yval)] <- ymax
xmax <- max(abs(xax))

# Define breaks in size legend based on data in the comparison
flog.info("Define point size") 
lowest <- xmax/4
medium <- max(xmax/2, lowest + 0.5) # Make sure to round to a higher number than lowest
highest <- xmax-0.5 # Make sure to round to lower number
size_breaks <- round(c(lowest, medium, xmax-0.5), digits = 0)
size_labels <- apply(expand.grid(size_breaks, as.vector(unlist(comparison[i]))), 1, paste, collapse=", ") # Set different values for different comparisons

flog.info("Define point shape")
# Variable used to set the shapes of the points (side of the plot)
sign_shape <- sign(xax)
sign_shape[sign_shape == 1] <- shapes[2]
sign_shape[sign_shape == -1] <- shapes[1]
sign_shape <- as.factor(sign_shape)

flog.info("Create and format labels") 
# Variable used to annotate genes in the plots
vollabels <- res$gene_name # Get gene names

# Add labels only in genes that are not hypothetical and have x > 3 and y > 80
vollabels[(abs(xax) < 3) & (yax < 80)] <- ""
vollabels[vollabels == "-"] <- ""

# Make labels bold and italic
italic_labels <- vollabels
italic_labels[!italic_labels == ""] <- paste0("bolditalic('", vollabels[!vollabels == ""],"')")
select_labs <- italic_labels[!italic_labels == ""]

flog.info("Assign colors to data points") 
# Variable to store the colors
keyvals.col <- c()
# Coloring depending on x (different colors if it's < 0.5, > 3 or between) and y
# (different colors depending on if the padj is <1e-5 or >1e-5)
keyvals.col <- ifelse(
  abs(res$log2FoldChange) < lfc_filter & res$padj < pval_filter, "#00008B",
  ifelse(abs(res$log2FoldChange) > lfc_filter & abs(res$log2FoldChange) < 3 & res$padj < pval_filter, "#008B8B",
         ifelse(abs(res$log2FoldChange) > 3 & res$padj < pval_filter, "#66CDAA",
                ifelse(abs(res$log2FoldChange) < lfc_filter & res$padj > pval_filter, "#8B7D6B",
                       ifelse(abs(res$log2FoldChange) > 3 & res$padj > pval_filter, "#EED5B7",
                              "#CDB79E")))))

# Labels for each color
keyvals.col[is.na(keyvals.col)] <- "red" # NAs should not be in the data
names(keyvals.col)[keyvals.col == "#66CDAA"] <- expression(italic("p"["adj"])*" < 0.1, Log"[2]*italic("FC")*" > 3")
names(keyvals.col)[keyvals.col == "#008B8B"] <- expression(italic("p"["adj"])*" < 0.1, 1 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "#00008B"] <- expression(italic("p"["adj"])*" < 0.1, Log"[2]*italic("FC")*" < 1")
names(keyvals.col)[keyvals.col == "#EED5B7"] <- expression(italic("p"["adj"])*" > 0.1, Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "#CDB79E"] <- expression(italic("p"["adj"])*" > 0.1, 1 < Log"[2]*italic("FC")*" < 3")
names(keyvals.col)[keyvals.col == "#8B7D6B"] <- expression(italic("p"["adj"])*" > 0.1, Log"[2]*italic("FC")*" < 1")

# Assign colors to data points
color_values <- unique(keyvals.col)
color_values <- color_values[order(color_values)] # Colors must be ordered to match the plot
labels <- unlist(lapply(color_values, function(color) {
  unique(names(keyvals.col)[keyvals.col == color])}))

# Vector to scale the size of data points
size_vector <- abs(res$log2FoldChange)

flog.info(paste("Generate the volcano plot for file"), input_files[i]) 
# Code to generate the Volcano plot
volcanoplot <- ggplot(data = res, aes(x = log2FoldChange, y = yval, col = keyvals.col, label = vollabels)) +
  geom_vline(xintercept = c(-lfc_filter, lfc_filter), col = "gray", linetype = "dashed") + # Add dashed line to show log2FC < 0.5
  geom_hline(yintercept = -log10(pval_filter), col = "gray", linetype = "dashed") + # Add dashed line for p-value > 1e-5
  geom_point(aes(size = size_vector, shape = sign_shape, fill = keyvals.col), alpha = 0.6, stroke = 0.5, color = "darkorchid") + # Line the size, shape, color (fill + border) and stroke of the points
  geom_label_repel(max.overlaps = Inf, color = "black", #show_guide = FALSE,  
        size = 3, box.padding = 0.4, fontface = "bold.italic") + # Add label boxes
  scale_shape_manual(values = shapes[order(shapes)], guide = "none") + # Set point shapes
  theme_light() + # Set point size and overall graph appearance
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # Remove grid
  scale_size_continuous(breaks = rep(size_breaks, times = 2), labels = size_labels, range = scale_val) + # Set the legend for point size
  scale_fill_manual(values = color_values, # Set the colors of up/downregulated points
      labels = parse(text = labels)) + # Set the color labels
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3, shape = 21, color = "darkorchid")), # Increase the size of legend points, change shape to add border color
       size = guide_legend(nrow = 3, ncol = 2, bycol = TRUE, order = 2,
       override.aes = list(shape = rep(shapes, each = 3), color = "darkorchid"))) + # Make size legend include shape information
  labs(fill = "Differential expression", # Set legend labels
       size = expression("Log"[2]*italic("Fold Change")), # Legend title
       x = expression("Log"[2]*italic("Fold Change")), # Title of main axes
       y = expression("-Log"[10]*italic("p"["adj"])), color = FALSE) +
  coord_cartesian(ylim = c(0, round(ymax)), xlim = c(-(round(xmax) + 1), round(xmax) + 1)) + # Set axis limits
  scale_x_continuous(breaks = seq(-(round(xmax) + round(xmax) + 1), 8, 1)) + # Customize ticks in the x axis
  scale_y_continuous(breaks = seq(0, round(ymax), 25)) + # Customize ticks in y axis
  ggtitle(title, subtitle = subtitle) + # Plot title
  theme(plot.title = element_text(hjust = 0.5), # Center title
       plot.subtitle = element_text(hjust = 0.5)) # Center subtitle

# Show plot in console
#volcanoplot

# Save the plots to files
flog.info(paste("Saving plot to ", output_files[2*i -1], " and to ", output_files[2*i]))
ggsave(output_files[2*i -1], width = 9, height = 6)
ggsave(output_files[2*i], width = 9, height = 6)
}

