# cfDNAmethy
Disease Prediction by Cell-Free DNA Methylation 

This is a basic example of using Quadratic Programming (QP) or Non-negative Matrix Factorization (NMF) on cfDNA methylation data. 

"analysis.R" is the main program. It requires "functions.R" and "dat.rda" to run. 

"functions.R" contains utility functions that need to be loaded in R before running "analysis.R". 
"dat.rda" contains two datasets:

(1) methylation level for all samples and features. It is a matrix of feature (row) by sample (column).

(2) external reference panel for all tissues. It is a matrix of feature (row) by tissue (column).

In the main program:

"build.QP" will use QP to estimate tissue proportions. It requires inputs of (1) methylation level for all samples and (2) tissue reference panel. 

"build.NMF" will use NMF to estimate both tissue reference panel and tissue proportions. It requires inputs of (1) methylation level for all samples and (2) matrix rank in deconvolution. 

Additionally, if user wants to evaluate the prediction accuracy, "build.model.svm.1cv" will conduct leave-one-out cross-validation (LOOCV) using Support Vector Machine (SVM). 
