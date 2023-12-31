"""
There was an issue with Bayespsace not clustering certain spots that were known to have UMIs 
I ran this script to try to figure out if the input matrix was missing these specific positions
or if it was the results of BayesSpace excluding them. 
Conclusion: The initial input matrix was missing these positions. 
This was a result of the filtering thresholds in SockEye pipeline. 
"""

library(SingleCellExperiment)
library(Matrix)
library(SpatialExperiment)
library(BayesSpace)
library(ggplot2)
library(readr)

ID <- 'Br2720_post'
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br2720_post/230707_Br2720_post_all_p10/expression_matrices/bayes_space_files'
CM <- '230707_Br2720_post_all_p10_gene_expression_processed.csv'
column_data <- '230707_Br2720_post_all_p10_gene_expression_processed_barcodes_flipped.csv'

CountsMatrix <- read.csv(paste0(path,'/', CM), row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
colData <- read.csv(paste0(path,'/', column_data), stringsAsFactors=FALSE, check.names=FALSE)


spe <- SingleCellExperiment(assays=list(logcounts=CountsMatrix), colData=colData)

set.seed(100)
dec <- scran::modelGeneVar(spe)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
spe <- scater::runPCA(spe, subset_row = top)

set.seed(102)
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE,
                         n.PCs=15, n.HVGs=2000, log.normalize=FALSE)
q <- 2
d <- 15

set.seed(103)
spe <- spatialCluster(spe, q=q, d=d, platform='Visium',
                      nrep=5000, gamma=3, save.chain=TRUE)
#labels <- dplyr::recode(spe$spatial.cluster, 1, 2, 3, 5, 8, 4, 7, 6, 9)
labels <- dplyr::recode(spe$spatial.cluster, 1, 2)
#labels <- dplyr::recode(spe$spatial.cluster, 1, 2, 3, 5, 4, 7, 6)

#clusterPlot(spe, palette=NULL, size=0.05) +
#  scale_fill_viridis_d(option = "A", labels = 1:2) +
#  labs(title="Br8667_mid_lib2")

p <- clusterPlot(spe, label=labels, palette=NULL, size=0.05) +
  scale_fill_brewer(palette = "Set1", labels = 1:2) +  # Adjusted this line
  labs(title=paste0(ID, '_2_clusters')) +
  coord_flip() + 
  scale_y_reverse()
p




# Load necessary library
library(readr)

# Define the file path
file_path <- "/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br2720_post/230707_Br2720_post_all_p10/expression_matrices/bayes_space_files/230707_Br2720_post_all_p10_gene_expression_processed_barcodes_flipped.csv"

# Read the data
data <- read_csv(file_path)

# Initialize a 78x128 matrix with zeros
grid <- matrix(0, nrow = 78, ncol = 128)

# Mark the positions from the file
for (i in 1:nrow(data)) {
  row <- data$row[i]
  col <- data$col[i]
  if (row <= 78 && col <= 128) {
    grid[row, col] <- 1
  }
}

# Identify missing positions
missing_positions <- which(grid == 0, arr.ind = TRUE)

# Print missing positions
print(missing_positions)

# Convert the matrix to a data frame for ggplot
grid_df <- as.data.frame(as.table(grid))

# Rename columns for clarity
names(grid_df) <- c("Row", "Col", "Value")

# Create the plot
ggplot(grid_df, aes(x = Col, y = Row, fill = factor(Value))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "red", "1" = "green")) +
    labs(fill = "Spot Status", x = "Column", y = "Row") +
    theme_minimal() +
    coord_fixed()
