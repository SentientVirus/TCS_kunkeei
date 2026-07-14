# Redirect all output to log file
log <- "~/mucoid_project/proteomics/logs/05a-volcano_proteomics.log"
con <- file(log, "w+")
sink(con, append = TRUE, type="message")
sink(con, append = TRUE)

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#

# Note: This is the only script for proteomic data analysis which requires
# the R environment
library("ggplot2")
library("ggrepel")
library(futile.logger)

#=============================================================================#
# 0. Logging                                                                  #
#=============================================================================#

# Create a logger that will be saved to a file
flog.logger("volcano", TRACE, appender=appender.file(log))

flog.info("R script to generate Volcano plots")

#=============================================================================#
# 1. Load dataframe with DESeq2 output + annotations                          #
#=============================================================================#

# Define inputs
flog.info("Definining input variables")
indir <- "~/mucoid_project/proteomics/Perseus/results/parsed"
outdir <- "~/mucoid_project/proteomics/plots"

input_files <- list.files(path = indir, pattern = ".tsv", full.names = TRUE)
output_files <- gsub("Perseus/results/parsed", "plots", input_files)
output_files <- gsub(".tsv", ".pdf", output_files)

#Function to change the first letter of the labels to uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Define vectors with formatting
titles <- c() # Plot titles
conditions <- c()
all_shapes <- c() # Desired point shape for the plot
comparison <- c() # Comparisons to be plotted
for (file in input_files){
  basefile <- basename(file)
  basefile <- substr(basefile, 1, nchar(basefile)-4)
  samples <- strsplit(basefile, "_")
  samples <- samples[[1]]
  cond1_1 <- substr(samples[1], 1, 3)
  cond1_2 <- substr(samples[1], 4, 4)
  cond1_3 <- samples[2]
  if (grepl("D", samples[1], fixed = TRUE)){
    cond1_4 <- "+D"
  }else{cond1_4 <- "-D"}
  cond2_1 <- substr(samples[4], 1, 3)
  cond2_2 <- substr(samples[4], 4, 4)
  cond2_3 <- samples[5]
  if (grepl("D", samples[4], fixed = TRUE)){
    cond2_4 <- "+D"
  }else{cond2_4 <- "-D"}
  titles <- c(titles, basefile)
  comparison <- c(comparison, basefile)
  if (cond1_1 == "muc" & cond2_1 == "agg"){
    shapes <- c(23, 21)
    conds <- c("Muc", "Agg")
  }
  else if (cond1_1 == "agg" & cond2_1 == "muc"){
    shapes <- c(21, 23)
    conds <- c("Agg", "Muc")
  }
  else if (cond1_2 == "F" & cond2_2 == "S"){
    shapes <- c(25, 24)
    conds <- c("-S", "+S")
  }
  else if (cond1_2 == "S" & cond2_2 == "F"){
    shapes <- c(24, 25)
    conds <- c("+S", "-S")
  }
  else if (cond1_3 == "log" & cond2_3 == "stat"){
    shapes <- c(21, 22)
    conds <- c("log", "stat")
  }
  else if (cond1_3 == "stat" & cond2_3 == "log"){
    shapes <- c(22, 21)
    conds <- c("stat", "log")
  }
  else if (cond1_4 == "+D" & cond2_4 == "-D"){
    shapes <- c(24, 25)
    conds <- c("+D", "-D")
  }
  else if (cond1_4 == "-D" & cond2_4 == "+D"){
    shapes <- c(25, 24)
    conds <- c("-D", "+D")
  }
  all_shapes[[length(all_shapes)+1]] <- shapes
  conditions[[length(conditions)+1]] <- conds
}
#shape_names <- rep(c("Morphology", "Sucrose"), each = 2)
scale_val <- c(0, 2.5) # Variable to scale points

# Loop through i to get res and plot
for (i in 1:length(input_files)){
  
# Read the annotated DESeq2 output
flog.info(paste("Reading file ", input_files[i])) 
res <- read.table(input_files[i], sep = "\t", 
                 numerals = "no.loss", header = TRUE, row.names = 1, quote = "",
                 stringsAsFactors = FALSE)

# Convert numeric columns to numeric (character by default)
j <- 7:17
res[, j] <- apply(res[, j], 2, function(x) as.numeric(x))

# Get shapes and titles for the plot
shapes <- as.vector(unlist(all_shapes[i]))
title <- titles[i]

# Get the mean log-scaled LFQs for each condition
cond1_avg <- rowMeans(res[, 12:14])
cond2_avg <- rowMeans(res[, 15:17])

#Remove the data points where both are 0 to avoid problems plotting 
#(point shapes are dependent on the condition where a gene is most expressed)
res <- res[cond1_avg != cond2_avg, ]

#=============================================================================#
# 2. Volcano plot with ggplot and standard plotting                           #
#=============================================================================#

# Set the x and y variables 
flog.info("Calculate x and y") 
yax <- res$Welch.s.T.test.q.value # Get the y values (q-values)
yax <- -log10(yax) # Calculate the -log10
cond1_avg <- rowMeans(res[, 12:14]) # Average log2 LFQ for the first condition
cond2_avg <- rowMeans(res[, 15:17]) # Average log2 LFQ for the second condition
xax <- cond1_avg-cond2_avg # Calculate the log2FC by subtracting the two log2 values

# Get maximum and set infinite values (padj = 0) to maximum
yval <- yax
ymax <- max(yval[is.finite(yval)])
yval[!is.finite(yval)] <- ymax
xmax <- max(abs(xax))

# Define breaks in size legend based on data in the comparison
flog.info("Define point size") 
lowest <- min(abs(xax))+0.1
medium <- max(xmax/2-0.1, lowest) # Make sure to round to a higher number than lowest
highest <- xmax-0.1 # Make sure to round to lower number

# Define the sizes to show in the legend
size_breaks <- round(c(lowest, medium, highest), digits = 1)
size_labels <- apply(expand.grid(size_breaks, as.vector(unlist(conditions[i]))), 1, paste, collapse=", ") # Set different values for different comparisons

flog.info("Define point shape")
# Variable used to set the shapes of the points (side of the plot)
sign_shape <- sign(xax)
sign_shape[sign_shape == 1] <- shapes[2]
sign_shape[sign_shape == -1] <- shapes[1]
sign_shape <- as.factor(sign_shape)

flog.info("Create and format labels") 
# Variable used to annotate genes in the plots
vollabels <- res$Refined.annotations # Get gene names
# Add labels only in genes that are not hypothetical and have x > 2 and y > 1
vollabels[(abs(xax) < 2) | (yax < -log10(0.05))] <- ""
vollabels[vollabels == "-"] <- ""

flog.info("Assign colors to data points") 
# Variable to store the colors
keyvals.col <- c()
# Coloring depending on x (different colors if it's < 0.5, > 3 or between) and y
# (different colors depending on if the padj is <1e-5 or >1e-5)
keyvals.col <- ifelse(
  abs(xax) < 2 & yval > -log10(0.05), "#00008B",
  ifelse(abs(xax) > 2 & abs(xax) < 4 & yval > -log10(0.05), "#008B8B",
         ifelse(abs(xax) > 4 & yval > -log10(0.05), "#66CDAA",
                ifelse(abs(xax) < 2 & yval < -log10(0.05), "#8B7D6B",
                       ifelse(abs(xax) > 4 & yval < -log10(0.05), "#EED5B7",
                              "#CDB79E")))))

# Labels for each color
keyvals.col[is.na(keyvals.col)] <- "red" # NAs should not be in the data
names(keyvals.col)[keyvals.col == "#66CDAA"] <- expression(italic("p"["adj"])*" < 0.05, Log"[2]*italic("FC")*" > 4")
names(keyvals.col)[keyvals.col == "#008B8B"] <- expression(italic("p"["adj"])*" < 0.05, 2 < Log"[2]*italic("FC")*" < 4")
names(keyvals.col)[keyvals.col == "#00008B"] <- expression(italic("p"["adj"])*" < 0.05, Log"[2]*italic("FC")*" < 2")
names(keyvals.col)[keyvals.col == "#EED5B7"] <- expression(italic("p"["adj"])*" > 0.05, Log"[2]*italic("FC")*" < 4")
names(keyvals.col)[keyvals.col == "#CDB79E"] <- expression(italic("p"["adj"])*" > 0.05, 2 < Log"[2]*italic("FC")*" < 4")
names(keyvals.col)[keyvals.col == "#8B7D6B"] <- expression(italic("p"["adj"])*" > 0.05, Log"[2]*italic("FC")*" < 2")

# Assign colors to data points
color_values <- unique(keyvals.col)
color_values <- color_values[order(color_values)] # Colors must be ordered to match the plot
labels <- unlist(lapply(color_values, function(color) {
  unique(names(keyvals.col)[keyvals.col == color])}))

# Vector to scale the size of data points
size_vector <- abs(xax)
flog.info(paste("Generate the volcano plot for file", input_files[i]))
# Code to generate the Volcano plot
volcanoplot <- ggplot(data = res, aes(x = xax, y = yval, col = keyvals.col, label = firstup(vollabels))) +
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = "dashed") + # Add dashed line to show log2FC < 0.5
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") + # Add dashed line for p-value > 1
  geom_point(aes(size = size_vector, shape = sign_shape, fill = keyvals.col), alpha = 0.6, stroke = 0.5, color = "darkorchid") + # Apply the size, shape, color (fill + border) and stroke of the points
  geom_label_repel(box.padding = 0.4, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE, color = "black", alpha = 0.7,
                   size = 3, fontface = "bold", seed = 1, nudge_x = -0.4, nudge_y = 0.2) + # Add label boxes
  scale_shape_manual(values = shapes[order(shapes)], guide = "none") + # Set point shapes
  theme_light() + # Set point size and overall graph appearance
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # Remove grid
  scale_size_continuous(breaks = rep(size_breaks, times = 2), labels = size_labels, range = scale_val) + # Set the legend for point size
  scale_fill_manual(values = color_values, # Set the colors of up/downregulated points
                    labels = parse(text = labels)) + # Set the color labels
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3, shape = 21, color = "darkorchid")), # Increase the size of legend points, change shape to add border color
         size = guide_legend(nrow = 3, ncol = 2, order = 2,
                             override.aes = list(shape = rep(shapes, each = 3), color = "darkorchid"))) + # Make size legend include shape information
  labs(fill = "Differential expression", # Set legend labels
       size = expression("Log"[2]*italic("Fold Change")), # Legend title
       x = expression("Log"[2]*italic("Fold Change")), # Title of main axes
       y = expression("-Log"[10]*italic("p"["adj"]))) +
  coord_cartesian(ylim = c(0, round(ymax+0.5)), xlim = c(-(round(xmax) + 1), round(xmax) + 1)) + # Set axis limits
  scale_x_continuous(breaks = seq(-(round(xmax) + round(xmax) + 1), 8, 1)) + # Customize ticks in the x axis
  scale_y_continuous(breaks = seq(0, round(ymax+0.5), 0.5)) + # Customize ticks in y axis
  ggtitle(title) + # Plot title
  theme(plot.title = element_text(hjust = 0.5), # Center title
        plot.subtitle = element_text(hjust = 0.5)) # Center subtitle

# Show plot in console
volcanoplot

# Save the plots to files
flog.info(paste("Saving plot to", output_files[i]))
ggsave(output_files[i], width = 9, height = 6)
}

