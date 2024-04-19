#=============================================================================#
# 0. Install packages if needed                                               #
#=============================================================================#
list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org');

to_install <- c("DESeq2", "ggplot2", "stringr", "apeglm", "pheatmap", "ggrepel")
new.packages <- to_install[!(to_install %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package);
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("stringr")
library(futile.logger)

# Script to generate a Volcano plot using ggplot, following this tutorial:
# https://biostatsquid.com/volcano-plots-r-tutorial/

yval <- -log10(res$padj)
yval[!is.finite(yval)] <- 300

sign_color <- sign(res$log2FoldChange)
sign_color[sign_color > 0] <- 0
sign_color[abs(res$log2FoldChange) < 0.5] <- 1
sign_color[res$padj > 0.05] <- 2
sign_color <- as.factor(sign_color)
vollabels <- rownames(res)
vollabels[(abs(res$log2FoldChange) < 4) & (-log10(res$padj) < 100)] <- ""

volcanoplot <- ggplot(data = res, aes(x = res$log2FoldChange, y = yval, col = sign_color, label = vollabels)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = "dashed") + # Add dashed line to show log2FC < 0.5
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") + # Add dashed line for p-value > 0.05
  geom_point(size = 2, alpha = 0.5, stroke = 0) + theme_light() + # Set point size and overall graph appearance
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # Remove grid
  scale_color_manual(values = c("#00AFBB", "#bb0c00", "black", "grey"), # Set the colors of up/downregulated points
                     labels = c("Inhibitor", "Mucoid", expression("Log"[2]*"FC < 0.5"), expression("p"["adj"]*" < 0.05"))) + # Set the color labels
  coord_cartesian(ylim = c(0, 320), xlim = c(-8, 8)) + # Set plot limits
  labs(color = "Differential expression", # Legend title
       x = expression("Log"[2]*"(Fold Change)"), y = expression("-Log"[10]*"(adjusted p-value)")) +
  scale_x_continuous(breaks = seq(-8, 8, 1)) + # Customize ticks in the x axis
  scale_y_continuous(breaks = seq(0, 320, 20)) + # Customize ticks in y axis
  ggtitle("Mucoid vs Inhibitor, sucrose") + # Plot title
  geom_text_repel(max.overlaps = Inf) # To show all labels

volcanoplot <- volcanoplot + 
  theme(plot.title = element_text(hjust = 0.5))
# Saving the plots to files
ggsave("try.png", width = 9, height = 6)