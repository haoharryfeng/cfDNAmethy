source("functions.R")
load("dat.rda")


########################
## use QP to estimate tissue proportions
########################
## input is (1) methylation level for all samples and (2) tissue reference panel
res1 = build.QP(mix, ref)
# res1 is the solved tissue proportion matrix, each column is a sample, each row is a tissue

########################
## use NMF to estimate both tissue reference panel and tissue proportions
########################
## input is (1) methylation level for all samples and (2) matrix rank in deconvolution
res2 = build.NMF(mix, rank = 14)
# res2 is a list object
# res2$W is the predicted W matrix, res2$H is the predicted H matrix


########################
## use leave-one-out cross-validation (LOOCV) to evaluate the prediction accuracy 
########################
#true.class is the ground truth for 
true.class = colnames(mix)

acc.marker = build.model.svm.1cv(true.class, t(mix)) #using marker directly
acc.QP = build.model.svm.1cv(true.class, t(res1)) #using QP
acc.NMF = build.model.svm.1cv(true.class, t(res2$H)) #using NMF



