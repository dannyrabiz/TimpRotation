library(SingleCellExperiment)
library(Matrix)
library(SpatialExperiment)
library(BayesSpace)
library(ggplot2)
library(readr)


control_data <- read.csv("control_file_transcript.csv")
for(i in 1:nrow(control_data)) {
  ID <- control_data[i, 'ID']
  path <- control_data[i, 'path']
  CM <- control_data[i, 'CountMat']
  column_data <- control_data[i, 'ColData']
  
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
                        nrep=50000, gamma=3, save.chain=TRUE)
  #labels <- dplyr::recode(spe$spatial.cluster, 1, 2, 3, 5, 8, 4, 7, 6, 9)
  labels <- dplyr::recode(spe$spatial.cluster, 1, 2)
  #labels <- dplyr::recode(spe$spatial.cluster, 1, 2, 3, 5, 4, 7, 6)
  
  #clusterPlot(spe, palette=NULL, size=0.05) +
  #  scale_fill_viridis_d(option = "A", labels = 1:2) +
  #  labs(title="Br8667_mid_lib2")
  
  p <- clusterPlot(spe, label=labels, palette=NULL, size=0.05) +
    scale_fill_brewer(palette = "Set1", labels = 1:2) +  # Adjusted this line
    labs(title=paste0(ID, 'transcripts_2_clusters')) +
    coord_flip() + 
    scale_y_reverse()
  p
  
  ggsave(paste0('2Clusters/Transcripts/', ID, '_k2.png'), plot = p, width = 10, height = 8, units = "in")
  
  # Assuming `dlpfc` is your SingleCellExperiment object and `spatial.cluster` is the column with cluster assignments
  cluster_assignments <- spe$spatial.cluster
  
  # Creating a data frame with barcodes and their corresponding cluster assignments
  output_data <- data.frame(Barcode = colnames(spe), Cluster = cluster_assignments)
  
  # Writing the data frame to a CSV file
  write.csv(output_data, paste0('2Clusters/Transcripts/',ID, '_k2_cluster_assignments.csv'), row.names = FALSE)
  }

