"""
Alternative pipeline for BayesSpace
Pipeline also referred to as 'Hope's Pipeline' in some of the comparison figures
Refer to the following file for more detail of each step
"""

library(SingleCellExperiment)
library(Matrix)
library(SpatialExperiment)
library(BayesSpace)
library(ggplot2)
library(readr)

CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667_mid_lib2/lib2_pass10/lib2_pass10.transcript_expression.processed.csv", row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667_mid_lib2/lib2_pass10/Br8667_lib2_transcript_processed_barcodes_coordinates_flipped.csv", row.names=1, header = TRUE)

spe <- SingleCellExperiment(assays=list(logcounts=CountsMatrix), colData=colData)

set.seed(100)
dec <- scran::modelGeneVar(spe)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
spe <- scater::runPCA(spe, subset_row = top)

set.seed(102)
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE,
                             n.PCs=15, n.HVGs=2000, log.normalize=FALSE)
q <- 9
d <- 15

set.seed(103)
spe <- spatialCluster(spe, q=q, d=d, platform='Visium',
                        nrep=5000, gamma=3, save.chain=TRUE)
labels <- dplyr::recode(spe$spatial.cluster, 1, 2, 3, 5, 8, 4, 7, 6, 9)


clusterPlot(spe, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A", labels = 1:9) +
  labs(title="Br8667_mid_lib2")

p <- clusterPlot(spe, label=labels, palette=NULL, size=0.05) +
  scale_fill_brewer(palette = "Set1", labels = 1:9) +  # Adjusted this line
  labs(title="9Clusters-Gene") +
  coord_flip() + 
  scale_y_reverse()
p

# same data with only 7 layers - like example from Bayes space vignette 
​
set.seed(104)
spe7 <- spatialCluster(spe, q=7, d=d, platform='Visium',
                        nrep=50000, gamma=3, save.chain=TRUE)
​
labels7 <- dplyr::recode(spe7$spatial.cluster, 3, 4, 5, 6, 2, 7, 1,)
​
​
clusterPlot(spe7, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A", labels = 1:7) +
  labs(title="BayesSpace")
