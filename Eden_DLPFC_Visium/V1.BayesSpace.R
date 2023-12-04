"""
This is the first pipeline I've used for running BayesSpace. 
Bayesspace's vignette tutorials have more detail on what is going on in each step
This pipeline does not allow for specification of # of clusters as easily as V2 - not exactly sure why yet. 
"""
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Matrix)

#import expression data
counts <- read.csv('/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667_mid_lib2/lib2_pass10/lib2_pass10.gene_expression.counts.csv', row.names=1, stringsAsFactors=FALSE, check.names=FALSE)

#import barcode data
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667_mid_lib2/lib2_pass10/Br8667_lib2_gene_processed_barcodes_coordinates_flipped.csv", stringsAsFactors=FALSE, check.names=FALSE)

counts_matrix <- as.matrix(counts)
counts_sparse <- as(counts_matrix, "dgCMatrix”)

#fix differences in barcounts and expression data
mismatch_in_counts <- setdiff(colnames(counts_sparse), colData$Barcode)
mismatch_in_colData <- setdiff(colData$Barcode, colnames(counts_sparse))

# Remove columns from counts_sparse that are not in colData
counts_sparse <- counts_sparse[, !(colnames(counts_sparse) %in% mismatch_in_counts)]
colData <- colData[!(colData$Barcode %in% mismatch_in_colData), ]

sce <- SingleCellExperiment(assays = list(counts = counts_sparse), colData = colData )

dlpfc <- spatialPreprocess(sce, platform="ST", n.PCs=7, n.HVGs=2000, log.normalize=TRUE)


dec <- scran::modelGeneVar(dlpfc)
top <- scran::getTopHVGs(dec, n = 2000)

q <- 7  # Number of clusters
d <- 15
#unclear if this step is necessary
dlpfc <- qTune(dlpfc, qs=seq(2, 10), platform="ST", d=7)

#50,000 is recommended - should take about 15 minutes
dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=3, save.chain=TRUE)


labels <- dplyr::recode(dlpfc$spatial.cluster, 3, 4, 5, 6, 2, 7, 1)

#Here is the flipped:
p <- clusterPlot(dlpfc, label=labels, palette=NULL, size=0.05) +
     scale_fill_brewer(palette = "Set1", labels = 1:7) +  # Adjusted this line
     labs(title="Br8667_mid_lib2") +
     coord_flip() + 
     scale_y_reverse()

#Save plot 
ggsave("Br8667_lib2.png", plot = p, width = 10, height = 8, units = "in")
