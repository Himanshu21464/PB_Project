library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesame)
library(SNFtool)
library(EnhancedVolcano)
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



# Perform SNF on the principal component matrices
snf_result <- SNFtool::SNF(list(gene_affinity, methylation_affinity), K, alpha)
view(snf_result)


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

cluster_assignments <- cluster_df$cluster


cluster_factor <- factor(cluster_assignments)



brca_matrix_com_df<-data.frame(brca_matrix_com)
colnames(cluster_rows_df)
library(limma)

library(edgeR)


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

view(results)
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

# Apply the function to each cluster
for (i in 1:7) {
  results <- calculate_logFC(i, results)
}




# Calculate the mean of the logFC columns
logFC_mean <- rowMeans(results[, grep("^logFC_", colnames(results))])

# Add a new column named "logfoldchange" to the results matrix
results$logfoldchange <- logFC_mean
View(results)
library(ggplot2)


# Subset the results matrix to include only columns of interest
results_sub <- results[, c("logfoldchange", "P.Value")]

# Set the threshold for significant genes
threshold <- exp(-60)

# Create a vector of colors for significant and non-significant genes
colors <- ifelse(results_sub$P.Value < threshold, "red", "black")

# Create the volcano plot
EnhancedVolcano(
  results_sub,
  lab = rownames(results_sub),
  x = "logfoldchange",
  y = "P.Value",
  xlim = c(-3, 3),
  title = "Enhanced Volcano Plot",
  pCutoff = threshold,
  FCcutoff = 0.0089,
  
)

# Filter the results to only include significant genes
significant_genes <- results_sub[results_sub$P.Value < threshold & results_sub$logfoldchange<0.0089, ]

# Add a new column with the gene names
significant_genes$GeneName <- rownames(significant_genes)
View(significant_genes)

view(significant_genes)
gene_list <- as.matrix(significant_genes[, c("GeneName", "logfoldchange")])
gene_list <- gene_list[order(gene_list[,2], decreasing=TRUE),]


#gene_list <- sort(gene_list, decreasing = TRUE)

# Extract the first 14 characters of the GeneName column
gene_list <- substr(gene_list, 1, 15)
library(ggrepel)
library(ggfortify)
library(plotly)

view(new_matrix)
view(gene_list)
library(dplyr)
library(stats)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
view(gene_list)
rownames(gene_list)


result <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   universe      = NULL)
print(result)


library(ggplot2)
library(RColorBrewer)
dotplot(result, showCategory = 30, title = "GO Enrichment Analysis" )
barplot(result, showCategory = 30, xlab= "GO Enrichment Analysis" )
