
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("maftools")
BiocManager::install("sesame")
BiocManager::install("proxy")
BiocManager::install("SNFtool")

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
# download and visualize mutation data from TCGA ----------------------
query_mutation <- GDCquery(project = 'TCGA-LAML',
                           data.category = 'Simple Nucleotide Variation',
                           access = 'open',
)

output_query_mutation <- getResults(query_mutation)

GDCdownload(query_mutation)

maf <- GDCprepare(query_mutation, summarizedExperiment = TRUE)
# mutation data is in maf
View(maf)
#maftools utils to read and create dashboard
maftools.input <- read.maf(maf)
View(maf)



plotmafSummary(maf = maftools.input,
               addStat = 'median',
               rmOutlier = TRUE,
               dashboard = TRUE)

#oncoprint
oncoplot(maf = maftools.input,
         top = 10,
         removeNonMutated = TRUE)

################################ Similarity Network Fusion ##########################################


# Mark the parameters
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)


# step 1 Data normalization
View(brca_matrix)
brca_matrix_scaled <- scale(brca_matrix)
brca_matrix_log <- log2(brca_matrix_scaled + 1)
gene_data <- na.omit(brca_matrix_log)
View(gene_data)


# Preprocess methylation data
methyl_data_scaled <- scale(methyl_data)
methyl_data_log <- log2(methyl_data_scaled + 1)
methyl_data_norm <- na.omit(methyl_data_log)

View(methyl_data_norm)

###########ERROR - 1 #############################
# Preprocess mutation data
# Error in normalize.quantiles(as.matrix(maf)) : 
# vector types do not match in copyVector
#maf_norm <- normalize.quantiles(as.matrix(maf))

###resolved###


maf <- maf[, colSums(is.na(maf)) == 0]
View(maf)
###########ERROR - 1 #############################


#gene_Dist1 = (dist2(gene_data,gene_data))^(1/2)
#methyl_Dist2 = (dist2(methyl_data_norm,methyl_data_norm))^(1/2)

## ERROR:
#> gene_Dist1 = (dist2(as.matrix(gene_data),as.matrix(gene_data)))^(1/2)
#Error: cannot allocate vector of size 27.4 Gb
#> gene_Dist1 = (dist2(gene_data,gene_data))^(1/2)
#Error: cannot allocate vector of size 27.4 Gb

## calculating euclidean/correlation distances


gene_data_dist <- proxy::dist(t(gene_data), method = "correlation")
#On viewing this it shows 
#> View(gene_data_dist)
#Error in names[[i]] : subscript out of bounds
View(gene_data_dist)

print(gene_data_dist)
# should be square matrix
print(dim(gene_data_dist))

#151 151
gene_data_dist_mat <- as.matrix(gene_data_dist)
dim(gene_data_dist_mat) <- c(151, 151)
colnames(gene_data_dist_mat) <-colnames(gene_data)
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
dim(methylation_dist_mat) <- c(194, 194)
colnames(methylation_dist_mat) <-colnames(methyl_data_norm)
class(methylation_dist_mat)
dim(methylation_dist_mat)
View(methylation_dist_mat)


## Calculating Affinity matrixes

gene_affinity = SNFtool::affinityMatrix(gene_data_dist_mat, K, alpha)
View(gene_affinity)
methylation_affinity = SNFtool::affinityMatrix(methylation_dist_mat, K, alpha)
View(methylation_affinity)


#performing SNF

#similarity_matrix = SNF(list(gene_affinity,methylation_affinity), K, T, merge.method = "intersection")

library(dplyr)

row.names(gene_data) <- 1:nrow(gene_data)
row.names(methyl_data_norm) <- 1:nrow(methyl_data_norm)

# Merge the two matrices based on their column names
merged_data <- inner_join(as.data.frame(gene_affinity), as.data.frame(methylation_affinity), by = "row.names")

# Remove the first column containing row names
merged_data <- merged_data %>% select(-1)

# Convert the merged data frame to a matrix
merged_matrix <- as.matrix(merged_data)

# Calculate the similarity matrix using SNF
similarity_matrix <- SNF(list(merged_matrix), K, T)












pheatmap(gene_data_dist)


# 151 151

View(methyl_data_log)


# Increase distance between columns

# 194 194

mutation_dist <- proxy::dist(t(maf), method = "correlation")
print(mutation_dist)
# 40 40
print(is.na(gene_data_dist)== FALSE)

## create Affinity matrixes

##########################ERROR-2##########################
#> gene_W1 = affinityMatrix(gene_data_dist, K, alpha)
#Error in x[cbind(i, i)] <- value : subscript out of bounds

#gene_data_sim = affinityMatrix(gene_data_dist, 20, 0.5)

# gene_data_sim <- SNFtool::affinityMatrix(gene_data_dist, K = 20, alpha = 0.5)
# methylation_sim <- SNFtool::affinityMatrix(methylation_dist, K = 20, alpha = 0.5)
# mutation_sim <- SNFtool::affinityMatrix(mutation_dist, K = 20, alpha = 0.5)

##########################ERROR-2##########################


##########################ERROR-3##########################

# similarity_gene <- computeSimilarity(t(brca_matrix_log))

# Error in computeSimilarity(t(brca_matrix_log)) : 
# could not find function "computeSimilarity"
##########################ERROR-3##########################


##########################ERROR-4##################################################
# gene_sim_mat <- proxy::simil(1/(1+brca_matrix_log+ 1e-8))
# Error in do.call(".External", c(list(CFUN, x, y, pairwise, if (!is.function(method)) get(method) else method),  : 
#                negative length vectors are not allowed
                                
##########################ERROR-4##################################################


##########################ERROR-5##################################################

# gene_data_sim <- SNFtool::affinityMatrix(gene_data_dist, K = 20)

# class(gene_data_dist)  -- "dist"
# > gene_data_sim <- SNFtool::SNF(list(gene_data_dist), K = 20)
# Error in rowSums(X) : 'x' must be an array of at least two dimensions


# > gene_data_sim <- affinityMatrix(gene_data_dist, K = 20, sigma = 0.5)
# Error in x[cbind(i, i)] <- value : subscript out of bounds
# gene_data_sim <- affinityMatrix(gene_data_dist, K = 20, sigma = 0.5)


##########################ERROR-5##################################################

