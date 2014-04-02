setwd("~/Desktop/Database/Rdata")
load("motstab.Rdata")

s1<-matrix(c(-1,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(-1,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)
s3<-matrix(c(-1,1,-1,-1,0,1,1,-1,0),nrow=3,ncol=3)
s4<-matrix(c(-1,0,1,0,0,1,-1,-1,0),nrow=3,ncol=3)
s5<-matrix(c(-1,1,1,-1,0,0,-1,0,0),nrow=3,ncol=3)
d1<-matrix(c(-1,1,1,1,0,1,0,0,0),nrow=3,ncol=3)
d2<-matrix(c(-1,1,1,0,0,1,0,1,0),nrow=3,ncol=3)
d3<-matrix(c(-1,1,1,1,0,0,0,0,0),nrow=3,ncol=3)
d4<-matrix(c(-1,1,1,0,0,0,0,1,0),nrow=3,ncol=3)
d5<-matrix(c(-1,1,1,0,0,1,1,0,0),nrow=3,ncol=3)
d6<-matrix(c(-1,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(-1,1,1,1,0,1,1,0,0),nrow=3,ncol=3)
d8<-matrix(c(-1,1,1,1,0,0,1,0,0),nrow=3,ncol=3)
s.bi<-matrix(c(-1,1,1,0,-1,0,0,1,-1,0,0,1,0,-1,-1,0),nrow=4,ncol=4)

mot.lst <- list(s1, s2, s3, s4, s5, d1, d2, d3, d4, d5, d6, d7, d8)
analyze.eigen <- function(m){
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- runif(1, 0, 10)
      }
      if(m[i, j] == -1){
        m[i, j] <- runif(1, -1, 0)
      }
    }
  }
  diag(m) <- -1
  ev <- max(Re(eigen(m)$values))
  return(ev)
}

test_stability <- function(matrices, rep){
  eigens<-lapply(matrices, function(x){replicate(rep, analyze.eigen(x))})
  qss <- lapply(eigens, function(x){sum(x < 0) / rep})
  return(qss)
}

system.time(
mot.stability <- test_stability(mot.lst, 1000000)
)
names(mot.stability) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
ms.mat <- unlist(mot.stability)
colnames(ms.mat) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
sort(ms.mat)


#-------------------------------------------------
get_qss<- function(chains, mode, parms, iter){
  # Input: list of matrices, uniform or normal dist, parameters, and number of iterations
  # Output: matrix of QSS values 
  # parameters is dataframe
  ## unif: min and max of distribution
  ## norm: mean and standard deviation
  test <- matrix(nrow = length(chains), ncol = nrow(parms))
  for(i in 1:nrow(parms)){
    eigen.test <- lapply(chains, analyze_eigen, mode = mode, iter = iter, 
                         params = parms[i,], self = 0)
    qss.test <- lapply(eigen.test, function(x){
      sum(x < 0) / 10000
    })
    test[,i] <- unlist(qss.test)
  }
  return(test)
}

eigen_lnorm <- function(m, params, self = -1){
  # For when I want to use normal distribution
  # Params is dataframe of mean and standard deviation for relative impact of prey on pred 
  ## and pred on prey
  ev <- c()
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- abs(rnorm(1, params$pred1, params$pred2))
        m[j, i] <- -abs(rnorm(1, params$prey1, params$prey2))
      }
    }
  }
  diag(m) <- self
  ev <- max(Re(eigen(m)$values))
  return(ev)
}

function(m, iter, mode, params, self = 0){
  # Input: matrix, number of iterations, unif or norm, parameters
  # Output: iter # of eigenvalues 
  if (mode == "unif"){
    evals <- c()
    for (i in 1:iter){
      eig <- eigen_unif(m, params, self = self)
      evals[i] <- eig
    }
    return(evals)
  }
  if (mode == "norm"){
    evals <- c()
    for (i in 1:iter){
      eig <- eigen_norm(m, params, self = self)
      evals[i] <- eig
    }
    if (mode == "lnorm"){
      evals <- c()
      for (i in 1:iter){
        eig <- eigen_lnorm(m, params, self = self)
        evals[i] <- eig
      }
    return(evals)
  }
}

library(devtools)
source_url("https://raw.github.com/jjborrelli/Ecological-Networks/master/Chain%20Length/Rscripts/chains_functions.R")

params.n <- data.frame(pred1 = c(0, 0, 0, 0, 0, 0, 0), pred2 = c(1, 2, 3, 4, 5, 6, 7), 
                       prey1 = c(0, 0, 0, 0, 0, 0, 0), prey2 = c(1, 1, 1, 1, 1, 1, 1))

mot.qss <- get_qss(mot.lst, mode = "norm", parms = params.n, iter = 10000)

params.u <- data.frame(pred1 = c(0, 0, 0, 0, 0, 0, 0), pred2 = c(10, 10, 10, 10, 5, 3, 1), 
                       prey1 = c(-5, -1, -.1, -.01, -1, -1, -1), prey2 = c(0, 0, 0, 0, 0, 0, 0))

mot.qss.u <- get_qss(mot.lst, mode = "unif", parms = params.u, iter = 10000)

