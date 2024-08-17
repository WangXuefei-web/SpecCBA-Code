### Functions for competitive methods ###
#' @import     rMultiNet # for HOOI

### Method: SpecCov ###
SpecCov.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecCov <- SpecCov(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecCov$loglik,alpha,re.SpecCov$NMI,re.SpecCov$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecCov$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)

  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'],
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

### Method: SpecCBA ###
SpecCBA.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecCBA <- SpecCBA(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecCBA$loglik,alpha,re.SpecCBA$NMI,re.SpecCBA$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecCBA$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'], 
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

SpecCBA.fast.fun <- function(sample,alpha.seq,K){
  start <- Sys.time()
  result.seq <- array(NA, dim=c(length(alpha.seq),4)) # loglik alpha NMI HamAcc
  lb.list <- list()
  alpha.i <- 1
  for (alpha in alpha.seq){
    re.SpecCBA <- SpecCBA(sample, alpha, K)
    result.seq[alpha.i,] <- c(re.SpecCBA$loglik,alpha,re.SpecCBA$NMI,re.SpecCBA$HamAcc)
    lb.list[[alpha.i]] <- list('alpha'=alpha,'label.hat'=re.SpecCBA$label.hat)
    alpha.i <- alpha.i+1
  }
  result.seq <- as.data.frame(result.seq)
  colnames(result.seq) <- c('loglik','alpha','NMI','HamAcc')
  max.id <- which.max(result.seq$loglik)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("alpha"=result.seq[max.id,'alpha'],"NMI"=result.seq[max.id,'NMI'], 
                 "HamAcc"=result.seq[max.id,'HamAcc'],"label.hat"=lb.list[[max.id]]$label.hat,"time"=as.numeric(time))
  return(output)
}

### Method: mean adjacency matrix ###
computeLaplacian <- function(Adj){
  d <- apply(Adj,1,sum)
  d[d==0] <- 1e-10
  D <- d^(0.5)
  D.inv <- D^(-1)
  G <- diag(D)
  G.inv <- diag(D.inv)
  L <- G.inv%*%Adj%*%G.inv
  return(L)
}

eigenratio_test <- function(mat){
  k_max <- floor((dim(mat)[1])/3)
  k_min <- min(3,dim(mat)[1])
  
  # spectral decomposition
  eigcov <- eigen(mat) 
  eigen_values <- eigcov$values # has sorted in decreasing order
  
  # eigenvalue ratio test for the number of factors
  eigen_ratios <- eigen_values[-length(eigen_values)] / eigen_values[-1]  
  eigen_ratios <- eigen_ratios[1:k_max]
  max_ratio_index <- max(k_min,which.max(eigen_ratios))
  
  # estimating loading matrix
  eigen_vectors <- eigcov$vectors[, 1:max_ratio_index]
  return(eigen_vectors)
}

MeanAdj.fun<-function(sample,K)
{
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # compute mean adjacency matrix and its Laplacian
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  Lap <- computeLaplacian(A.bar)

  eigenvectors <- eigenratio_test(Lap)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(eigenvectors, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  
  return(output)
}

### Method: OLMF orthogonal linked matrix factorization ###
## The objective function 

lmffunctiono<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  objloop<- sum(unlist(lapply(1:M,function(m){
    specobj<-norm(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar),type="F")^2
    return(specobj)
  })))
  obj=objloop
  return(obj)
}


##  The gradients

lmfdero<-function(param,laplist,n,k){
  M=length(laplist)
  ustar<-matrix(param[1:(n*k)],n,k)
  lambda<-lapply(1:M,function(m){return(matrix(param[(n*k+(m-1)*k^2+1):(n*k+m*k^2)],k,k))})
  derlist1<-lapply(1:M,function(m){
    specobj= -(diag(n)-ustar%*%t(ustar))%*%laplist[[m]]%*%ustar%*%lambda[[m]]
    return(specobj)
  })
  derlist2<-lapply(1:M,function(m){
    specobj= -t(ustar)%*%(laplist[[m]]-ustar%*%lambda[[m]]%*%t(ustar))%*%ustar
    return(specobj)
  })
  der1<-Reduce("+",derlist1)
  der2<-unlist(derlist2)
  return(c(as.vector(der1),as.vector(der2)))
}


## The main function with BFGS optimization

OLMF.fun<-function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # compute the Laplacian of each layer
  Laplist <- lapply(1:L,function(l){return(computeLaplacian(A@data[,,l]))})
  
  # Initialize with mean Laplacian
  Lap.bar <- Reduce('+',Laplist)
  spectra <- eigen(Lap.bar)
  ustar <- spectra$vectors[,1:K]
  lambda <- lapply(1:L,function(l){return(diag(spectra$values[1:K]))})
  
  # Optimization
  param <- c(as.vector(ustar),as.vector(unlist(lambda)))
  optimized <- optim(par=param,fn=lmffunctiono,gr=lmfdero,method="BFGS",control=list(reltol=0.0001,maxit=200),laplist=Laplist,n=N,k=K)
  param <- optimized$par
  
  # Community detection
  ustar <- matrix(param[1:(N*K)],N,K)
  lambda <- lapply(1:L,function(l){return(matrix(param[(N*K+(l-1)*K^2+1):(N*K+l*K^2)],K,K))})
  
  specstar <- kmeans(ustar,K)
  label.hat <- specstar$cluster
  
  # Compute metrics
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- end - start
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  
  return(output)
  
}

### Method: SOSBA ###
SOSBA.fun <- function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  # computations based on adjacency tensor
  A.l <- as.array(A@data) # N*N*L
  
  # debias
  A.bar <- as.matrix(modeMean(A,3,drop=TRUE)@data)
  d.bar <- colSums(A.bar)
  D.bar <- diag(d.bar)
  
  # compute covariance-type matrix
  AA <- meanXXT(A.l,transpose=FALSE)
  S <- L * (AA - D.bar)
  
  # estimating U
  U <- loading_mat_by_eigenratio(S,K)
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(U, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  # metric
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}




### Method: HOOI higher order orthogonal iteration ###
HOOI.fun <- function(sample,K){
  start <- Sys.time()
  A <- sample$adjT
  N <- dim(A)[1]
  L <- dim(A)[3]
  label <- sample$global_membership
  
  tucker.re <- tucker(A,ranks=c(K,K,1))
  U <- tucker.re$U[[1]]
  
  # community detection 
  get.est.flag <- FALSE
  temp <- list()
  try.temp <- 1
  while(!get.est.flag){
    temp <- try(kmeans(U, K, nstart = 10),silent=TRUE)
    if('try-error' %in% class(temp)) # judge weather error occurs
    {
      try.temp <- try.temp + 1
      if (try.temp == 3)
      {
        stop("There seems to be an issue with the initialization within the kmeans step. Please run the algorithm again.") 
      }
      next
    }else{
      get.est.flag <- TRUE
    }
  }
  kmeans.re <- temp 
  label.hat <- kmeans.re$cluster
  
  # metric
  NMI <- nmi(label,label.hat)
  HamAcc <- MCrate(label, label.hat)
  
  end <- Sys.time()
  time <- difftime(end,start,units = "secs")
  output <- list("NMI"=NMI, "HamAcc"=HamAcc,"time"=as.numeric(time))
  return(output)
}

