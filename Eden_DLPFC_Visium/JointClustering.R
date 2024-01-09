library(BayesSpace)

#Part 1 
#Import Data

#UNFILTERED SAMPLES:

#Sample1 - Br8667
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667mid_all/br8667mid_p10/expression_matrices/bayes_space_files'
CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667mid_all/br8667mid_p10/expression_matrices/bayes_space_files/br8667mid_p10.gene_expression_counts.csv", row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667mid_all/br8667mid_p10/expression_matrices/bayes_space_files/br8667mid_p10_gene_barcodes_coordinates_flipped.csv", row.names=1, header = TRUE)
spe0 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#Sample2 - Bra 720
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br2720_post/230707_Br2720_post_all_p10/expression_matrices/bayes_space_files'
CM <- '230707_Br2720_post_all_p10_gene_expression_counts.csv'
column_data <- '230707_Br2720_post_all_p10_gene_expression_counts_barcodes_flipped.csv'

CountsMatrix <- read.csv(paste0(path,'/', CM), row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
colData <- read.csv(paste0(path,'/', column_data), stringsAsFactors=FALSE, check.names=FALSE)
spe1 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#Sample3 - Br6522
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files'
CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files/230707_Br6522_ant_all_p10.gene_expression.counts.csv"
                         , row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files/230707_Br6522_ant_all_p10.gene_expression.counts_barcodes_flipped.csv",
                    row.names=1, header = TRUE)
spe2 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#Sample4 - Br6471
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files'
CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files/230707_Br6471_post_all_p10.gene_expression.counts.csv"
                         , row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files/230707_Br6471_post_all_p10.gene_expression.counts_barcodes_flipped.csv",
                    row.names=1, header = TRUE)
spe3 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#------------------------------------------------------------------------------
#Processed/FILTERED

#Sample1

CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667mid_all/br8667mid_p10/expression_matrices/bayes_space_files/br8667mid_p10.gene_expression_processed.csv", row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_br8667_mid/br8667mid_all/br8667mid_p10/expression_matrices/bayes_space_files/br8667mid_p10_gene_processed_barcodes_coordinates_flipped.csv", row.names=1, header = TRUE)
spe0 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)


#Sample2 - Bra 2720
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br2720_post/230707_Br2720_post_all_p10/expression_matrices/bayes_space_files'
CM <- '230707_Br2720_post_all_p10_gene_expression_processed.csv'
column_data <- '230707_Br2720_post_all_p10_gene_expression_processed_barcodes_flipped.csv'

CountsMatrix <- read.csv(paste0(path,'/', CM), row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
colData <- read.csv(paste0(path,'/', column_data), stringsAsFactors=FALSE, check.names=FALSE)
spe1 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#Sample3 - Br6522
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files'
CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files/230707_Br6522_ant_all_p10.gene_expression.processed.csv"
                         , row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6522_ant/230707_Br6522_ant_all_p10/expression_matrices/bayes_space_files/230707_Br6522_ant_all_p10.gene_expression.processed_barcodes_flipped.csv",
                    row.names=1, header = TRUE)
spe2 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)


#Sample4 - Br6471
path <- '/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files'
CountsMatrix <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files/230707_Br6471_post_all_p10.gene_expression.processed.csv"
                         , row.names=1, check.names=F, stringsAsFactors=FALSE)
colData <- read.csv("/wilee/Data/eden/libd/sockeye_analysis/DLPFC_Br6471_post/230707_Br6471_post_all_p10/expression_matrices/bayes_space_files/230707_Br6471_post_all_p10.gene_expression.processed_barcodes_flipped.csv",
                    row.names=1, header = TRUE)
spe3 <- SingleCellExperiment(assays=list(counts=CountsMatrix), colData=colData)

#---------------------------------------------------------------------
#Part 2
#Format and Join single cell objects

# List of gene vectors from each SingleCellExperiment object
gene_lists <- list(rownames(spe0), rownames(spe1), rownames(spe2), rownames(spe3))



# Find common genes
common_genes <- Reduce(intersect, gene_lists)


# Subset both SCE objects to only include common genes
spe0 <- spe0[common_genes, ]
spe1 <- spe1[common_genes, ]
spe2 <- spe2[common_genes, ]
spe3 <- spe3[common_genes, ]

# Make sure the genes are in the same order
spe0 <- spe0[order(rownames(spe0)), ]
spe1 <- spe1[order(rownames(spe1)), ]
spe2 <- spe2[order(rownames(spe2)), ]
spe3 <- spe3[order(rownames(spe3)), ]


#Add missing 'Barcode' column if it doesn't exist
if(!"Barcode" %in% colnames(colData(spe0))) {
  colData(spe0)$Barcode <- NA
}

if(!"Barcode" %in% colnames(colData(spe1))) {
  colData(spe1)$Barcode <- NA
}

if(!"Barcode" %in% colnames(colData(spe2))) {
  colData(spe2)$Barcode <- NA
}

if(!"Barcode" %in% colnames(colData(spe3))) {
  colData(spe3)$Barcode <- NA
}


colData(spe0)$sample_name <- "Br8667mid"
colData(spe1)$sample_name <- "Br2720_post"
colData(spe2)$sample_name <- "Br6522_ant"
colData(spe3)$sample_name <- "Br6471_post"

#Combine into 1 SCE and preprocess
sce.combined = cbind(spe0, spe1, spe2,spe3, deparse.level = 1)


if (!"counts" %in% assayNames(sce.combined)) {
  names(assays(sce.combined))[1] <- "counts"
}

#sce.combined <- spatialPreprocess(sce.combined, platform="Visium", skip.PCA=FALSE,  n.PCs=15, n.HVGs=2000, log.normalize=TRUE)
#sce.combined = spatialPreprocess(sce.combined,  n.PCs = 15)#n.PCs = 50, platform="Visium", skip.PCA=TRUE, n.HVGs=2000, log.normalize=FALSE) #lognormalize, PCA
#sce.combined = spatialPreprocess(sce.combined, n.PCs = 50) 
sce.combined <- spatialPreprocess(sce.combined, platform="Visium", skip.PCA=FALSE,  n.PCs=15, n.HVGs=2000, log.normalize=TRUE)

###-------------------------------------------
#Part 3
#BATCH CORRECTION

sce.combined = runUMAP(sce.combined, dimred = "PCA")
colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")

#Original
ggplot(data.frame(reducedDim(sce.combined, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

ggplot(data.frame(reducedDim(sce.combined, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, 
           color = factor(sce.combined$sample_name), 
           shape = factor(sce.combined$sample_name))) +  # Mapping both color and shape
  geom_point(alpha = 0.1) +  # Adjust alpha for transparency
  labs(color = "Sample", shape = "Sample") +  # Label for color and shape legends
  theme_bw()


#RUN HARMONY
sce.combined = RunHarmony(sce.combined, "sample_name", verbose = F)
sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
  geom_point(alpha = 0.5) +
  labs(color = "Sample") +
  theme_bw()


#With Harmony
sce.combined = spatialCluster(sce.combined, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY

#Without Harmony 
sce.combined = spatialCluster(sce.combined, q=7, d=50, platform='Visium',
               nrep=5000, gamma=3, save.chain=TRUE)



sce_sample0 <- sce.combined[, sce.combined$sample_name == "Br8667mid"]
sce_sample1 <- sce.combined[, sce.combined$sample_name == "Br2720_post"]
sce_sample2 <- sce.combined[, sce.combined$sample_name == "Br6522_ant"]
sce_sample3 <- sce.combined[, sce.combined$sample_name == "Br6471_post"]


# Plot for Sample1
p1 <- clusterPlot(sce_sample0, color = NA) +
  labs(title = "Br8667mid") +
  theme_bw()


# Plot for Sample2
p2 <- clusterPlot(sce_sample1, color = NA) +
  labs(title = "Br2720_post") +
  theme_bw()

p3 <- clusterPlot(sce_sample2, color = NA) +
  labs(title = "Br6522_ant") +
  theme_bw()


p4 <- clusterPlot(sce_sample3, color = NA) +
  labs(title = "Br6471_post") +
  theme_bw()

#p1 + p2 + p3 + p4 


combined_plot <- p1 + p2 + p3 + p4 + 
  plot_layout(nrow = 2, ncol = 2)

ggsave("combined_plot.png", combined_plot, width = 10, height = 8, dpi = 300)


cluster_assignments_sample1 <- sce_sample1$spatial.cluster
output_data_sample1 <- data.frame(Barcode = colnames(sce_sample1), Cluster = cluster_assignments_sample1)
write.csv(output_data_sample1, "ClusterAssignment_Sample1.csv", row.names = FALSE)

cluster_assignments_sample2 <- sce_sample2$spatial.cluster
output_data_sample2 <- data.frame(Barcode = colnames(sce_sample2), Cluster = cluster_assignments_sample2)
write.csv(output_data_sample2, "ClusterAssignment_Sample2.csv", row.names = FALSE)



