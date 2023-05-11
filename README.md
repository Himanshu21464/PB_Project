# Practical Bioinformatics project for Signature Discovery Using Multi Omics Data Approach in Leukemia using unsupervised learning


# Code Readme

This repository contains R code for retrieving and analyzing gene expression and DNA methylation data from the TCGA-LAML project using the TCGAbiolinks package. The code performs various data retrieval, preprocessing, and analysis steps, including:

1. Retrieving project information and summary for TCGA-LAML.
2. Building a query to retrieve gene expression data.
3. Downloading gene expression data.
4. Preparing gene expression data for analysis.
5. Building a query to retrieve DNA methylation data.
6. Downloading DNA methylation data.
7. Preprocessing and analyzing DNA methylation data.
8. Performing Similarity Network Fusion (SNF) analysis on gene expression and DNA methylation data.
9. Clustering and visualizing the results of SNF analysis.
10. Differential expression analysis using limma.
11. Volcano plot visualization of differentially expressed genes.
12. GO enrichment analysis of differentially expressed genes.
13. Visualization of GO enrichment analysis results using dotplot, barplot, cnetplot, and goplot.
14. Creating a gene-gene network based on GO enrichment analysis results.

## Dependencies

The code requires the following R packages to be installed:

- TCGAbiolinks
- tidyverse
- maftools
- pheatmap
- SummarizedExperiment
- sesame
- SNFtool
- EnhancedVolcano
- limma
- ggrepel
- ggfortify
- plotly
- edgeR
- preprocessCore
- proxy
- dplyr
- stats
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- RColorBrewer
- enrichplot
- igraph

You can install these packages using the `install.packages()` function in R.

## Usage

1. Install the required R packages mentioned above.
2. Run the code in an R environment or an R script.
3. The code will retrieve gene expression and DNA methylation data from the TCGA-LAML project, preprocess the data, perform SNF analysis, and generate various visualizations and analysis results.
4. The results include gene expression data, DNA methylation data, clustering results, differential expression analysis results, volcano plots, GO enrichment analysis results, and gene-gene network visualization.

Please note that this code is specifically designed for the TCGA-LAML project and may need modifications to work with other TCGA projects or datasets.


## Steps

Library Imports: The code begins by importing various R packages that are required for data retrieval, preprocessing, analysis, and visualization. These packages include TCGAbiolinks, tidyverse, maftools, pheatmap, SummarizedExperiment, sesame, SNFtool, EnhancedVolcano, limma, ggrepel, ggfortify, plotly, edgeR, preprocessCore, proxy, dplyr, stats, clusterProfiler, org.Hs.eg.db, ggplot2, RColorBrewer, enrichplot, and igraph. These packages provide functions and tools for working with TCGA data, data manipulation, visualization, statistical analysis, and more.

Get Project Summary: The code uses the TCGAbiolinks package to retrieve information and summary of the TCGA-LAML project. The getGDCprojects() function retrieves a list of available TCGA projects, and the getProjectSummary() function retrieves the summary information for the specified project (TCGA-LAML in this case).

Build Query for Gene Expression Data: The code builds a query using the GDCquery() function to retrieve gene expression data from the TCGA-LAML project. The query specifies the project, data category (Transcriptome Profiling), experimental strategy (RNA-Seq), workflow type (STAR - Counts), and access type (open).

Retrieve Gene Expression Data: The getResults() function is used to retrieve the gene expression data based on the query constructed in the previous step.

Download Gene Expression Data: The GDCdownload() function is used to download the gene expression data files to the local machine.

Prepare Gene Expression Data: The code uses the GDCprepare() function to preprocess and prepare the gene expression data for analysis. The function creates a SummarizedExperiment object, which stores the gene expression data along with relevant metadata.

Build Query for DNA Methylation Data: Similar to step 3, this part of the code builds a query to retrieve DNA methylation data from the TCGA-LAML project. The query specifies the project, data category (DNA Methylation), platform (Illumina Human Methylation 27), and access type (open).

Retrieve DNA Methylation Data: The getResults() function is used to retrieve the DNA methylation data based on the query constructed in the previous step.

Download DNA Methylation Data: The GDCdownload() function is used to download the DNA methylation data files to the local machine.

Preprocess and Analyze DNA Methylation Data: The code performs preprocessing steps on the DNA methylation data using the GDCprepare() function. The resulting object contains the preprocessed methylation data.

Similarity Network Fusion (SNF): SNF is a technique used to integrate multiple similarity matrices to create a consensus similarity matrix. In this code, SNF is used to integrate the gene expression data and DNA methylation data. The integration is performed in the following steps:

First, similarity matrices are calculated for both the gene expression data and DNA methylation data using the affinityMatrix() function from the SNFtool package. The similarity matrices represent the pairwise similarity between samples based on the gene expression and DNA methylation profiles.

The affinityMatrix() function takes as input the distance matrices (gene_data_dist_mat and methylation_dist_mat) calculated from the gene expression and DNA methylation data, respectively, along with the number of nearest neighbors (K) and a hyperparameter (alpha).

The result is two affinity matrices: gene_affinity and methylation_affinity, which capture the similarity between samples based on gene expression and DNA methylation, respectively.

Clustering: Once the affinity matrices are obtained, clustering is performed to group samples based on their integrated similarities using the spectral clustering algorithm. The steps involved are:

The spectralClustering() function from the SNFtool package is used for spectral clustering. It takes the integrated similarity matrix (snf_result) as input and the number of desired clusters (C).

The spectralClustering() function assigns each sample to a specific cluster, resulting in a cluster assignment vector.

The displayClusters() function from the SNFtool package is used to visualize the clustering results, where each cluster is represented by a different color.

The cluster assignment vector (cluster) is converted to a data frame (cluster_df) for further processing and analysis.

Cluster Analysis: The cluster analysis involves organizing the samples and their corresponding gene expressions based on the identified clusters. The steps are:

An empty list (cluster_rows) is created to store the samples' row names corresponding to each cluster.

Using a loop, each cluster is iterated, and the indices of samples belonging to the current cluster are extracted from the cluster_df.

The row names of the samples in the current cluster are obtained from the snf_result_df and stored in the cluster_rows list.

After the loop, cluster_rows contains a list of row names for each cluster, where each element of the list represents a cluster and contains the corresponding sample names.

The results are stored in the cluster_rows_df data frame, which has two columns: cluster representing the cluster number and gene representing the sample names.


Differential Expression Analysis: The code uses the limma package, which provides functions for linear modeling and differential expression analysis, to perform differential expression analysis on the gene expression data. The lmFit() function is used to fit a linear model to the gene expression data, and the eBayes() function is used to perform empirical Bayes moderation of the log-fold changes. Finally, the topTable() function is used to extract the results, including log-fold changes and p-values, for all genes.

Volcano Plot Visualization: The code uses the EnhancedVolcano package to create a volcano plot, which is a commonly used visualization for differential expression analysis. The EnhancedVolcano() function is called, and it takes the results from the differential expression analysis as input. The volcano plot visualizes the log-fold changes on the x-axis and the negative log10 p-values on the y-axis. Genes that are significantly differentially expressed are highlighted with red color.

GO Enrichment Analysis: The code performs Gene Ontology (GO) enrichment analysis on the differentially expressed genes. GO enrichment analysis identifies the overrepresentation of GO terms, which represent biological processes, molecular functions, and cellular components, among the differentially expressed genes. The enrichGO() function from the clusterProfiler package is used for this analysis. It takes the gene list (significant genes) as input, along with the relevant parameters such as the organism database (org.Hs.eg.db), key type (ENSEMBL), ontology (BP for biological process), p-value cutoff, and q-value cutoff. The result is a list of enriched GO terms.

GO Term Visualization: The code visualizes the GO enrichment results using various plot types provided by the enrichplot package. The dotplot(), barplot(), cnetplot(), and goplot() functions are used to generate different types of plots to represent the enriched GO terms. These plots provide insights into the biological processes that are significantly associated with the differentially expressed genes.

Network Visualization: The code creates a network visualization using the igraph package to show the relationship between the enriched GO terms and the corresponding genes. It constructs a graph from the edges data frame, where each edge represents a connection between a GO term and a gene. The plot() function is then used to visualize the network graph.

These steps collectively perform differential expression analysis, GO enrichment analysis, and visualization of the results, allowing for the identification of biologically meaningful processes and genes associated with the analyzed data.





## Acknowledgments

This code utilizes the TCGAbiolinks package and other related R packages for data retrieval, preprocessing, analysis, and visualization. The TCGAbiolinks package provides a convenient interface to access and analyze TCGA data. Please refer to the respective package documentation and publications for more information.

## Disclaimer

This code is provided as-is without any warranty. Use it at your own risk.


