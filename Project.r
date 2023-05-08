
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
pheatmap(methyl_data[idx,])

View(methyl_data)
         

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
pheatmap(gene_affinity)
pheatmap(methylation_affinity)





library(dplyr)
library(stats)


pca_gene <- prcomp(gene_affinity, scale. = TRUE)
pca_methylation <- prcomp(methylation_affinity, scale. = TRUE)
dim
# Extract the first few principal components
pc_gene <- pca_gene$x[, 1:2]
pc_methylation <- pca_methylation$x[, 1:2]

dim(pc_gene)
dim(pc_methylation)
View(methylation_affinity_reduced)
K <- 20
alpha <- 0.5

# Perform SNF on the principal component matrices
snf_result <- SNFtool::SNF(list(gene_affinity, methylation_affinity), K, alpha)
view(snf_result)
num_clusters<-5

clustering_result <- spectralClustering(snf_result, num_clusters)

# View the clustering result
print(clustering_result)

