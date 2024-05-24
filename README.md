This is a guid to use and execute the GSEA code.

Input:
1. 3 file paths to: pathway gene sets, gene expression dataset, phenotypes
2. C: class of distinction in the alphabetic order of the phenotypes (ex: 0 for female and 1 for male)
3. p: wight of step in the running sum (0 for a Kolmogorov–Smirnov statistic, 1 for using the exact value of the fold change)
4. plotAllpathways: False to compute the enrichment score and its p-value of a specifed pathway, True to plot a list of pathways ESs with their adjusted p-values (for multiple hypothesis)

Procedure:
1. The program ranks the gene expression dataset
2. if plotAllpathways is true: compute ES and p-value for a specified pathway in the variable "pathway"
3. otherwise, plot all the pathways specified in the variable pathways
