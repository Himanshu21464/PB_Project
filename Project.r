
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesame)
library(SNFtool)


library(preprocessCore)
library(proxy)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LAML')


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-LAML',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
)

getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)


# prepare data
tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'fpkm_unstrand')
# gene data: brca_matrix
View(brca_matrix)

# build a query to retrieve DNA methylation data --------------
query_methly <- GDCquery(project = 'TCGA-LAML',
                         data.category = 'DNA Methylation',
                         platform = 'Illumina Human Methylation 27',
                         access = 'open',
                         data.type = 'Methylation Beta Value',
                         #barcode = c('TCGA-AB-2941-03A-01D-0743-05', 'TCGA-AB-2838-03A-01D-0741-05'))
                         
)

output_query_methyl <- getResults(query_methly)

GDCdownload(query_methly)

dup_samples <- duplicated(colnames(query_methly))
query_methly <- query_methly[, !dup_samples]
query_methly

#View()
# plot probes showing differences in beta values between samples                                               
dna.meth <- GDCprepare(query_methly, summarizedExperiment = TRUE)
methyl_data = assay(dna.meth)  
View(methyl_data)
#methylation data is in methyl_data

idx <- dna.meth %>% 
  assay %>% 
  rowVars() %>% 
  order(decreasing = TRUE) %>% 
  head(50)


# plot
pheatmap(methyl_data[idx,],fontsize_row =3,fontsize_col = 3)

#View(methyl_data)
         

################################ Similarity Network Fusion ##########################################


# Mark the parameters
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

meth_pat=substr(colnames(methyl_data),1,15)
gene_pat=substr(colnames(brca_matrix),1,15)
common_pat=intersect(meth_pat,gene_pat)
colnames(brca_matrix)  = substr(colnames(brca_matrix),1,15)
colnames(methyl_data)  = substr(colnames(methyl_data),1,15)
meth_com = methyl_data[,common_pat]
brca_matrix_com = brca_matrix[,common_pat]
# step 1 Data normalization
View(brca_matrix)
brca_matrix_scaled <- scale(brca_matrix_com)
brca_matrix_log <- log2(brca_matrix_scaled + 1)
gene_data <- na.omit(brca_matrix_log)
View(gene_data)
## Abhi tak toh hai.

# Preprocess methylation data
methyl_data_scaled <- scale(meth_com)
methyl_data_log <- log2(methyl_data_scaled + 1)
methyl_data_norm <- na.omit(methyl_data_log)

View(methyl_data_norm)




##################calculating distance matrix######################
gene_data_dist <- proxy::dist(t(gene_data), method = "correlation")
#On viewing this it shows 
#> View(gene_data_dist)
#Error in names[[i]] : subscript out of bounds
View(gene_data_dist)

print(gene_data_dist)
# should be square matrix
print(dim(gene_data_dist))
dim(gene_data)
dim(gene_data_dist)

#151 151

gene_data_dist_mat <- as.matrix(gene_data_dist)
#dim(gene_data_dist_mat) <- c(151, 151)
dim(gene_data_dist_mat)
#colnames(gene_data_dist_mat) <-colnames(gene_data)
View(gene_data_dist_mat)

methylation_dist <- proxy::dist(t(methyl_data_norm), method = "Euclidean")
#On viewing this it shows 
#> View(methylation_dist)
#Error in names[[i]] : subscript out of bounds

View(methylation_dist)

print(methylation_dist)
# should be square matrix
print(dim(methylation_dist))
# 194 194
methylation_dist_mat<-as.matrix(methylation_dist)

class(methylation_dist_mat)
dim(methylation_dist_mat)
View(methylation_dist_mat)


## Calculating Affinity matrixes

gene_affinity = SNFtool::affinityMatrix(gene_data_dist_mat, K, alpha)
View(gene_affinity)
methylation_affinity = SNFtool::affinityMatrix(methylation_dist_mat, K, alpha)
View(methylation_affinity)

dim(gene_affinity)
dim(methylation_affinity)
pheatmap(gene_affinity,fontsize_col = 3,fontsize_row = 3)
pheatmap(methylation_affinity,fontsize_col = 3,fontsize_row = 3)


# Perform SNF on the principal component matrices
snf_result <- SNFtool::SNF(list(gene_affinity, methylation_affinity), K, alpha)
view(snf_result)

# Create a function to calculate the average silhouette width
calculate_silhouette_width <- function(data, num_clusters) {
  # Perform spectral clustering
  cluster_labels <- spectralClustering(data, num_clusters)
  # Calculate the silhouette width
  silhouette_width <- silhouette(cluster_labels, dist(data))
  silhouette_width <- as.data.frame(silhouette_width)
  return(mean(silhouette_width$sil_width))
}

# Define the range of the number of clusters to test
num_clusters_range <- 2:8

# Calculate the silhouette width for each number of clusters
silhouette_widths <- sapply(num_clusters_range, function(num_clusters) {
  calculate_silhouette_width(snf_result, num_clusters)
})

# Plot the silhouette widths against the number of clusters
plot(num_clusters_range, silhouette_widths, type = "b", 
     xlab = "Number of clusters", ylab = "Average silhouette width")

## 7 is highest avg width in given range
## do clustering for 7 

#C = 7
#cluster<- spectralClustering(snf_result,C)
#displayClusters(snf_result,cluster)

#cluster_df <- as.data.frame(cluster)
#snf_result_df <- as.data.frame(snf_result)

#c1 <- which(cluster_df$cluster==1, arr.ind=TRUE)
#c1_row <- rownames(snf_result_df[c1, ])

#c2 <- which(cluster_df$cluster==2, arr.ind=TRUE)
#c2_row <- rownames(snf_result_df[c2, ])

#c3 <- which(cluster_df$cluster==3, arr.ind=TRUE)
#c3_row <- rownames(snf_result_df[c3, ])

#c4 <- which(cluster_df$cluster==4, arr.ind=TRUE)
#c4_row <- rownames(snf_result_df[c4, ])

#c5 <- which(cluster_df$cluster==5, arr.ind=TRUE)
#c5_row <- rownames(snf_result_df[c5, ])

#c6 <- which(cluster_df$cluster==6, arr.ind=TRUE)
#c6_row <- rownames(snf_result_df[c6, ])

#c7 <- which(cluster_df$cluster==7, arr.ind=TRUE)
#c7_row <- rownames(snf_result_df[c7, ])

################### Differential Expression Analysis ######################
C<-7

cluster<- spectralClustering(snf_result,C)
displayClusters(snf_result,cluster)

cluster_df <- as.data.frame(cluster)
snf_result_df <- as.data.frame(snf_result)

# Create an empty list to store the results
cluster_rows <- list()

# Loop over each cluster
for (i in 1:7) {
  # Get the indices of the data points in the current cluster
  c <- which(cluster_df$cluster == i, arr.ind = TRUE)
  
  # Extract the row names of the data points in the current cluster
  c_row <- rownames(snf_result_df[c, ])
  
  # Store the row names in the list
  cluster_rows[[i]] <- c_row
}

# View the results
cluster_rows

cluster_rows_df <- data.frame(cluster = rep(1:length(cluster_rows), sapply(cluster_rows, length)),
                              gene = unlist(cluster_rows))

# Extract the cluster assignments from the cluster_df data frame
cluster_assignments <- cluster_df$cluster

# Convert cluster_assignments into a factor variable
cluster_factor <- factor(cluster_assignments)



brca_matrix_com_df<-data.frame(brca_matrix_com)
colnames(cluster_rows_df)
library(limma)

library(edgeR)

# Convert the gene expression matrix into a DGEList object
dge <- DGEList(counts = brca_matrix_com_df)


view(dge)


view(cluster_factor)

design <- model.matrix(~ 0 + cluster_factor)
view(design)

v <- voom(dge, design)
view(v)

fit <- lmFit(v, design)

fit <- eBayes(fit)
results <- topTable(fit, number = Inf)

# View the results
view(results)



########################  logFC calculation #######################################


# Define a function to calculate logFC for each cluster
calculate_logFC <- function(cluster_name, results) {
  # Extract the column corresponding to the given cluster
  cluster_col <- paste0("cluster_factor", cluster_name)
  cluster_vals <- results[, cluster_col]
  
  # Calculate the average expression of the cluster
  cluster_aveexpr <- results$AveExpr
  
  # Calculate the logFC
  logFC <- log2(cluster_vals / cluster_aveexpr)
  
  # Add the logFC column to the results matrix
  results[, paste0("logFC_", cluster_name)] <- logFC
  
  return(results)
}


# Calculate the mean of the logFC columns
logFC_mean <- rowMeans(results[, grep("^logFC_", colnames(results))])

# Add a new column named "logfoldchange" to the results matrix
results$logfoldchange <- logFC_mean
View(results)


library(ggplot2)

# create a subset of the results matrix with columns of interest
results_sub <- results[, c("logfoldchange", "P.Value")]

# define the threshold for significant genes
threshold <- 0.05

# create a volcano plot using ggplot2
volcano_plot <- ggplot(results_sub, aes(x = logfoldchange, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value < threshold, "sig", "not_sig")), alpha = 0.7) +
  scale_color_manual(values = c("sig" = "black", "not_sig" = "red"), guide = FALSE) +
  labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot") +
  theme_bw()

# display the plot
volcano_plot

###################################### Volcano plot seems to be incorrect ###########################################
