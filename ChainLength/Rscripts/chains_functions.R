list_chains <- function(start, length){
  # Generates a list of matrices for universally omnivorous food chains
  # start is the length of the shortest chain
  # length is the length of the longest chain
  # creates chains of lengths between start and length
  chains <- list()
  for(i in start:length){
    chains[[i]] <- matrix(0, nrow = i, ncol = i)
  }
  
  for(i in start:length){
    chains[[i]][which(upper.tri(chains[[i]]))] <- 1
  }
  return(chains[start:length])
}

eigen_unif <- function(m, params, self = -1){
  # For when I want to use uniform distribution
  # Params is dataframe of mean and standard deviation for relative impact of prey on pred 
  ## and pred on prey
  ev <- c()
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- runif(1, params$pred1, params$pred2)
        m[j, i] <- runif(1, params$prey1, params$prey2)
      }
    }
  }
  diag(m) <- self
  ev <- max(Re(eigen(m)$values))
  return(ev)
}


eigen_norm <- function(m, params, self = -1){
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

eigen_lnorm <- function(m, params, self = -1){
  # For when I want to use normal distribution
  # Params is dataframe of mean and standard deviation for relative impact of prey on pred 
  ## and pred on prey
  ev <- c()
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- abs(rlnorm(1, params$pred1, params$pred2))
        m[j, i] <- -abs(rlnorm(1, params$prey1, params$prey2))
      }
    }
  }
  diag(m) <- self
  ev <- max(Re(eigen(m)$values))
  return(ev)
}

analyze_eigen <- function(m, iter, mode, params, self = -1){
  # Input: matrix, number of iterations, unif or norm, parameters
  # Output: iter # of eigenvalues 
  if (mode == "unif"){
    evals <- c()
    for (i in 1:iter){
      eig <- eigen_unif(m, params, self = -1)
      evals[i] <- eig
    }
    return(evals)
  }
  if (mode == "norm"){
    evals <- c()
    for (i in 1:iter){
      eig <- eigen_norm(m, params, self = -1)
      evals[i] <- eig
    }
    return(evals)
  }
  if (mode == "lnorm"){
    evals <- c()
    for (i in 1:iter){
      eig <- eigen_lnorm(m, params, self = -1)
      evals[i] <- eig
    }
    return(evals)
  }
}


find_qss <- function(chains, mode, parms, iter){
  # Input: list of matrices, uniform or normal dist, parameters, and number of iterations
  # Output: matrix of QSS values 
  # parameters is dataframe
  ## unif: min and max of distribution
  ## norm: mean and standard deviation
  test <- matrix(nrow = length(chains), ncol = nrow(parms))
  for(i in 1:nrow(parms)){
    eigen.test <- lapply(chains, analyze_eigen, mode = mode, iter = iter, 
                           params = parms[i,])
    qss.test <- t(sapply(eigen.test, function(x){
      sum(x < 0) / 10000
    }))
    test[,i] <- qss.test
  }
  return(test)
}