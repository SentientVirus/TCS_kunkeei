list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

other.packages <- c("NOISeq") #, "metaSeq")
new.packages <- other.packages[!(other.packages %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package)
}
#IDEA: No filtering according to TPM, filter here

library(NOISeq)
#library(EDASeq)
directory <- "featureCounts_reverse/countfiles/filtered" #"../featureCounts_reverse/countfiles/filtered"
parent_dir <- "featureCounts_reverse/countfiles" #"../featureCounts_reverse/countfiles"
metadat <- grep("metadata.tsv", list.files(parent_dir), value = TRUE)
plot_dir <- "plots"

if (!file.exists(plot_dir)){
  dir.create(plot_dir)
}

sampleFiles <- list.files(directory)[grepl("counts.tsv",list.files(directory))][-1]
sampleCondition <- c()
for (file in sampleFiles){
  if ((grepl("01", file, fixed = TRUE) == 1) | (grepl("02", file, fixed = TRUE) == 1)){
    sampleCondition <- c(sampleCondition, "mucoid")
  }
  else{sampleCondition <- c(sampleCondition, "inhibitor")}
}

sampleTable <- data.frame(sampleName = substring(sampleFiles, first = 10, last = 11),
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

filename <- sprintf("%s/%s", directory, sampleFiles[1])
counts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
for (file in sampleFiles[2:length(sampleFiles)]){
  filename <- sprintf("%s/%s", directory, file)
  cts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid"))
  counts <- cbind(counts, cts)
}

# Get unique substrings
unique_substrings <- unique(substr(colnames(counts), 1, 5))

# Initialize an empty dataframe to store the sums
sum_counts <- data.frame(row.names = rownames(counts))

# Iterate through the unique substrings and calculate sums
for (substring in unique_substrings) {
  matching_columns <- colnames(counts)[grepl(substring, colnames(counts))]
  sums <- rowSums(counts[, matching_columns])
  sum_counts[substring] <-  sums
}

# Rename the columns in the sums_dataframe
colnames(sum_counts) <- unique_substrings

# Now back to the analysis
coldata <- read.csv(sprintf("%s/%s", parent_dir, metadat), sep = "\t", row.names = 1)

# Do the same to coldata
# First, remove the 'condition3' column
grouped_coldata <- coldata[, -which(names(coldata) == "condition3")]

grouped_coldata$substring <- substr(rownames(coldata), 1, 5)

# Create a new dataframe with one row per unique combination of 'condition1' and 'condition2'
grouped_coldata <- grouped_coldata[!duplicated(grouped_coldata[c("substring", "condition1", "condition2")]), ]

# Assign the "substring" column as the row names
rownames(grouped_coldata) <- grouped_coldata$substring

# Remove the "substring" column
grouped_coldata$substring <- NULL

# Read chromosome names as metadata
chr <- read.delim("meta/gene_positions.tsv", sep = "\t")
rownames(chr) <- chr$Geneid
chr$Geneid <- NULL

chr2 <- chr[rownames(chr) %in% rownames(sum_counts), ]
chr <- chr[rownames(chr) %in% rownames(counts), ]

mydata2 <- readData(data=sum_counts, chromosome=chr2, factors=grouped_coldata)
postscript(file="plots/saturation_collapsed_k0.ps")
sat <- dat(mydata2, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:8)
dev.off()
png(file="plots/saturation_collapsed_k0.png", width=1000, height=750)
sat <- dat(mydata2, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:8)
dev.off()

mydata <- readData(data=counts, chromosome=chr, factors=coldata)
postscript(file="plots/saturation_k0.ps")
sat <- dat(mydata, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:40)
dev.off()
png(file="plots/saturation_k0.png", width=1000, height=750)
sat <- dat(mydata, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:40)
dev.off()