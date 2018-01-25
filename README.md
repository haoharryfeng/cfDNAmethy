# cfDNAmethy
Disease Prediction by Cell-Free DNA Methylation 

This is a basic example of using Quadratic Programming (QP) or Non-negative Matrix Factorization (NMF) on cfDNA methylation data. 

"analysis.R" is the main program. It requires "functions.R" and "dat.rda" to run. 

"functions.R" contains utility functions that need to be loaded in R before running "analysis.R". 
"dat.rda" contains two datasets:

(1) methylation level for all samples and features. It is a matrix of feature (row) by sample (column).

(2) external reference panel for all tissues. It is a matrix of feature (row) by tissue (column).

