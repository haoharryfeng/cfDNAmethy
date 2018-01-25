##### some utility functions
library(ROCR)
library(quadprog)
library(e1071)


#using Quadratic Programming to estimate proportions
#Y is the sample matrix of N_marker*N_people, M is the reference matrix of N_marker*N_tissue
build.QP <- function(Y, M) {
  p_est = matrix(0, ncol=ncol(Y), nrow = ncol(M))
  for(i in 1:ncol(Y)){
    dvec = t(M) %*% Y[,i]
    Dmat = t(M) %*% M
    Amat = cbind(rep(1, ncol(M)),diag(1,ncol(M)))
    bvec = c(1,rep(0,ncol(M)))
    meq = 1
    
    res=solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    p_est[,i]=res$solution
  }
  return(p_est)
}


build.NMF <- function(Y, rank) {
  out = RefFreeCellMethy(Y/100, mu0=RefFreeCellMixInitialize(Y/100, K = rank))
  H.pred = t(out$Omega)
  W.pred = out$Mu
  output = list(W = W.pred, H = H.pred)
  return(output)
}

### Methy
RefFreeCellMethy <- function(Y, mu0=NULL, K=NULL,iters=50,Yfinal=NULL,verbose=TRUE){
  if(is.null(mu0)){
    if(K==1) {
      if(!is.null(Yfinal)) Y <- Yfinal
      n <- dim(Y)[2]
      
      mu <- matrix(apply(Y,1,mean,na.rm=TRUE),ncol=1)
      omega <- matrix(1, n, 1)
      o <- list(Mu=mu, Omega=omega)
      class(o) <- "RefFreeCellMix"
      return(o)
    }
    else mu0 <- RefFreeCellMixInitialize(Y,K=K,method="ward")
  }
  for(i in 1:iters){
    flag <- !apply(is.na(mu0),1,any)
    omega <- projectMethy(Y[flag,],mu0[flag,])           # d=w2^T%*%E
    mu <- projectMethy(t(Y), omega, sumLessThanOne=FALSE) # d=H%*%Y^T
    if(verbose) print(summary(abs(as.vector(mu-mu0))))
    mu0 <- mu
  }
  if(!is.null(Yfinal)){
    mu <- projectMethy(t(Yfinal),omega,sumLessThanOne=FALSE)
  }
  
  o <- list(Mu=mu, Omega=omega)
  class(o) <- "RefFreeCellMix"
  o
}


####
projectMethy <- function(Y, Xmat, nonnegative=TRUE, sumLessThanOne=TRUE, lessThanOne=!sumLessThanOne){
  
  nCol = dim(Xmat)[2]   # the number of k
  nSubj = dim(Y)[2]
  
  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)
  
  if(nonnegative){
    # H matrix
    if(sumLessThanOne){
      Amat = cbind(rep(1,nCol), rep(-1,nCol), diag(nCol))
      b0vec = c(1, -1,rep(0,nCol))
      
      # W, calculate H
      for(i in 1:nSubj){
        obs = which(!is.na(Y[,i]))
        Dmat = t(Xmat[obs,])%*%Xmat[obs,]
        mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
      }
    }
    # W matrix
    else if(lessThanOne){
      Amat = cbind(-diag(nCol), diag(nCol))
      b0vec = c(rep(-1,nCol),rep(0,nCol))
      
      # H, calculate W
      for(i in 1:nSubj){
        obs = which(!is.na(Y[,i]))
        Dmat = t(Xmat[obs,])%*%Xmat[obs,]
        mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
      }
    }
    else{
      Amat = diag(nCol)
      b0vec = rep(0,nCol)
    }
    
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }
  
  return(mixCoef)
}
# Y: matrix appearing in the NMF to be factorized
# Xmat:The W or H matrix after NMF matrix Factorization
# Dmat:matrix appearing in the quadratic function to be minimized.
# dvec:vector appearing in the quadratic function to be minimized.
# Amat:matrix defining the constraints under which we want to minimize the quadratic function.
# bvec:vector holding the values of b0 (defaults to zero)

RefFreeCellMixInitialize <- function(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL, 
                                     largeOK=FALSE, dist.method = "euclidean", ...){
  
  if(!is.matrix(Y) | !is.numeric(Y)){
    stop("Y is not a numeric matrix\n")
  }
  n <- dim(Y)[2]
  
  if(is.null(Y.Cluster)){
    if(is.null(Y.Distance)){
      if(n>2500 & !largeOK){
        stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
      }
      Y.Distance <- dist(t(Y),method=dist.method)
    }
    Y.Cluster <- hclust(Y.Distance,...)
  } 
  
  classes <- cutree(Y.Cluster, K)
  s <- split(1:n,classes)
  
  sapply(s, function(u) apply(Y[,u,drop=FALSE],1,mean,na.rm=TRUE))
}

###############################################
## leave one out cross validation
###############################################

build.model.svm.1cv <- function(Y, X) {
  acc<- rep(0, length(Y))
  predY <- rep(NA, length(Y))
  for(i in 1:length(Y)) {
    cat(i, ",")
    Y.train <- Y[-i];    X.train <- X[-i,,drop=FALSE]
    model <- svm(x=X.train, y=factor(Y.train), probability=TRUE)
    predY[i] <- predict(model, X[i,,drop=F])
  }
  tbl<- table(Y, predY)
  tmp <- sum(diag(tbl)) / sum(tbl)
  res <- max(tmp, 1-tmp)
  return(res)
}





