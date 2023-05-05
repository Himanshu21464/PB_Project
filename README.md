# Practical Bioinformatics project for Signature Discovery Using Multi Omics Data Approach in Leukemia using unsupervised learning


This code performs data analysis on the TCGA-LAML project, which is a part of the TCGA (The Cancer Genome Atlas) project that aims to improve our understanding of the molecular basis of cancer. 

### Libraries used

This code requires the following libraries to be installed and loaded:

```
BiocManager
TCGAbiolinks
tidyverse
maftools
pheatmap
SummarizedExperiment
sesame
SNFtool
preprocessCore
proxy
```

You can install these packages using the `BiocManager::install()` function.

### Data Retrieval and Preparation

The code retrieves gene expression, DNA methylation, and mutation data from TCGA-LAML using the `GDCquery` function of `TCGAbiolinks`. 

After retrieving the data, the `GDCdownload` function downloads it. The `GDCprepare` function is then used to prepare the data, which returns a `SummarizedExperiment` object. 

The gene expression data is saved in the `brca_matrix` variable, the DNA methylation data is saved in the `methyl_data` variable, and the mutation data is saved in the `maf` variable.

### Data Analysis

#### DNA Methylation Analysis

The code generates a heatmap to visualize the DNA methylation data using the `pheatmap` library. 

#### Mutation Analysis

The code generates a summary plot and an oncoprint for the mutation data using the `maftools` library.

#### Similarity Network Fusion (SNF)

The code performs a similarity network fusion analysis on the gene expression and DNA methylation data using the ‘SNF’ function from the `SNFtool` library. 

Before the analysis, the gene expression and DNA methylation data are preprocessed by scaling and log2 transforming the data. The `proxy::dist()` function calculates the distance matrix, and then we transformed the matrices into square matrices to calculate the affinity matrix’ for gene expression and DNA methylation data using the “SNFtool::affinityMatrix()” function.

After performing the above steps, we got two matrices, “gene_affinity” and “methylation_affinity,” with dimensions 151x151 and 194x194, respectively.

The input matrices must be of the exact dimensions to perform Similarity Network Fusion.
So to reduce the dimensions of the “methylation_affinity” matrix from 194x194 to 151x151, we had to perform Principle Component Analysis (A statistical technique that is used to reduce the dimensionality of a dataset while retaining as much of its variance as possible.)

Then the SNF algorithm is applied to the list of 'Affinity matrices' to obtain the fused similarity matrix. 


### Errors and Resolutions

The code throws some errors, which have been resolved and explained in the comments. 

### Parameters

The parameters used for SNF are:

```
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)
```


