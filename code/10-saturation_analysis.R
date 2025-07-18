# Redirect all output to log file
con <- file(snakemake@log[[1]], "a+")
sink(con, append = TRUE, type="message")
sink(con, append = TRUE)

#=============================================================================#
# 0. Install packages if needed                                               #
#=============================================================================#

list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

other.packages <- c("NOISeq")
new.packages <- other.packages[!(other.packages %in% installed.packages()[,"Package"])]
for (package in new.packages){
  BiocManager::install(package)
}

#=============================================================================#
# 0. Load required libraries                                                  #
#=============================================================================#

library(NOISeq)
library(futile.logger)

#=============================================================================#
# 0. Logging                                                                  #
#=============================================================================#

# Create a logger that will be saved to a file
flog.logger("saturation", TRACE, appender=appender.file(snakemake@log))

flog.info("R script to run a saturation analysis")

#=============================================================================#
# 1. Load input variables from Snakemake                                      #
#=============================================================================#

flog.info("Loading inputs...")
directory <- dirname(snakemake@input[["counts"]][1]) #Directory with gene counts
metadat <- snakemake@input[["meta"]] #File with metadata
plot_dir <- dirname(snakemake@output[["plots"]][1]) #Directory to store the plots
sampleFiles <- snakemake@input[["counts"]] #Files with gene counts

# Create output directories if they don't exist
if (!file.exists(plot_dir)){
  dir.create(plot_dir)
}

#=============================================================================#
# 2. Loop through count files to retrieve the counts                          #
#=============================================================================#

flog.info("Retrieving counts...")
filename <- sampleFiles[1] #Get the first input file
counts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid")) #Create a matrix with the counts
for (filename in sampleFiles[2:length(sampleFiles)]){ #Loop through the remaining files
  cts <- as.matrix(read.csv(filename, sep="\t", row.names="Geneid")) #Create a matrix with counts
  counts <- cbind(counts, cts) #Add the matrix to the previous count matrix
}

#=============================================================================#
# 3. Calculate the mean values of counts per isolate/condition                #
#=============================================================================#

flog.info("Collapsing replicates...")

# Retrieve unique isolate/condition headers
unique_substrings <- unique(substr(colnames(counts), 1, 5))

# Initialize an empty dataframe to store the sums
sum_counts <- data.frame(row.names = rownames(counts))

# Iterate through the unique substrings and calculate sums
for (substring in unique_substrings) {
  matching_columns <- colnames(counts)[grepl(substring, colnames(counts))]
  sums <- rowSums(counts[, matching_columns])
  sum_counts[substring] <-  sums
}

#=============================================================================#
# 4. Collapse metadata                                                        #
#=============================================================================#

flog.info("Collapsing metadata...")

# Rename the columns in the sums_dataframe
colnames(sum_counts) <- unique_substrings

# Read metadata
coldata <- read.csv(metadat, sep = "\t", row.names = 1)

# Remove the last column (sequencing date)
grouped_coldata <- coldata[, -which(names(coldata) == "condition3")]

# Create a variable to collapse the metadata
grouped_coldata$substring <- substr(rownames(coldata), 1, 5)

# Create a new dataframe with one row per unique isolate/condition
grouped_coldata <- grouped_coldata[!duplicated(grouped_coldata[c("substring", 
                                             "condition1", "condition2")]), ]

# Assign the "substring" column as the row names
rownames(grouped_coldata) <- grouped_coldata$substring

# Remove the "substring" column
grouped_coldata$substring <- NULL


#=============================================================================#
# 5. Read position information (not needed)                                   #
#=============================================================================#

flog.info("Reading positional information...")

# Read chromosome names as metadata
chr <- read.delim(snakemake@input[["pos"]], sep = "\t") #Read metadata file
rownames(chr) <- chr$Geneid #Use locus tags as row index
chr$Geneid <- NULL #Remove the column with locus tags

# Remove genes that have been filtered out
chr2 <- chr[rownames(chr) %in% rownames(sum_counts), ]
chr <- chr[rownames(chr) %in% rownames(counts), ]

#=============================================================================#
# 6. Generate saturation plots                                                #
#=============================================================================#

flog.info("Generating saturation plots...")

# Run the saturation analysis collapsing replicates
mydata2 <- readData(data=sum_counts, chromosome=chr2, factors=grouped_coldata) #Read sums of counts
postscript(file=snakemake@output[["ps"]][1]) #Create .ps plot
sat <- dat(mydata2, type = "saturation", factor = NULL, k = 0, ndepth = 20) #Apply the saturation analysis without furthering filtering the data, get 20 datapoints
explo.plot(sat, samples = 1:8) #Plot the 8 isolate/condition combinations as one line each
invisible(dev.off()) #Restart the canvas
png(file=snakemake@output[["plots"]][1], width=1000, height=750) #Do the same, but to generate a PNG output
sat <- dat(mydata2, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:8)
invisible(dev.off())

# Run the saturation analysis for each sample
mydata <- readData(data=counts, chromosome=chr, factors=coldata)
postscript(file=snakemake@output[["ps"]][2])
sat <- dat(mydata, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:40) #The 40 samples are plotted
invisible(dev.off())
png(file=snakemake@output[["plots"]][2], width=1000, height=750)
sat <- dat(mydata, type = "saturation", factor = NULL, k = 0, ndepth = 20)
explo.plot(sat, samples = 1:40)
invisible(dev.off())

flog.info("Analysis completed!")
