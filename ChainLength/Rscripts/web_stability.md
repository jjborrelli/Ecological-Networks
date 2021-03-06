Niche model food webs and stability
========================================================

Load required libraries

```r
require(devtools)
```

```
## Loading required package: devtools
```

```r
require(igraph)
```

```
## Loading required package: igraph
```

```r
require(NetIndices)
```

```
## Loading required package: NetIndices
## Loading required package: MASS
```


Use the niche model function created by Ted Hart of Distributed Ecology

```r
niche.model <- function(S, C) {
    new.mat <- matrix(0, nrow = S, ncol = S)
    ci <- vector()
    niche <- runif(S, 0, 1)
    r <- rbeta(S, 1, ((1/(2 * C)) - 1)) * niche
    
    for (i in 1:S) {
        ci[i] <- runif(1, r[i]/2, niche[i])
    }
    
    r[which(niche == min(niche))] <- 1e-08
    
    for (i in 1:S) {
        
        for (j in 1:S) {
            if (niche[j] > (ci[i] - (0.5 * r[i])) && niche[j] < (ci[i] + 0.5 * 
                r[i])) {
                new.mat[j, i] <- 1
            }
        }
    }
    
    new.mat <- new.mat[, order(apply(new.mat, 2, sum))]
    return(new.mat)
}
```


Set the desired levels of number of nodes and connectances

```r
nodes <- seq(10, 100, 10)
connectance <- seq(0.05, 0.45, 0.025)
```


Create a data frame of node/connectance pairings

```r
niche.inputs <- data.frame(N = rep(nodes, each = 17), C = rep(connectance, 10))
```


Generate 10 niche model food webs for each pairing of node and connectance

```r
web_maker <- function(n, x) {
    # where n is # of webs to make for each pairing and x is a data frame of
    # node/connectance pairings Make sure the column names in the data frame are
    # N and C
    if (!colnames(x[1]) == "N") {
        warning("Check column names")
    }
    if (!colnames(x[2]) == "C") {
        warning("Check column names")
    }
    
    web <- list()
    
    for (i in 1:length(x[, 1])) {
        subweb <- list()
        for (j in 1:n) {
            subweb[[j]] <- niche.model(x$N[i], x$C[i])
        }
        web[[i]] <- subweb
    }
    
    return(web)
}
```


Create a list of niche model food webs for each pair of parameters

```r
system.time(test1 <- web_maker(10, niche.inputs))
```

```
##    user  system elapsed 
##  29.689   0.098  29.941
```


Call eigenvalue stability analysis functions from GitHub repository

```r
source("~/Desktop/GitHub/Ecological-Networks/ChainLength/Rscripts/chains_functions.R")
```


Define the parameters for `find_qss`

```r
params.n <- data.frame(pred1 = 0, pred2 = 2, prey1 = 0, prey2 = 1)
```


Analyze the quasi-sign stability of the niche model food webs with `find_qss`

```r
system.time(test2 <- lapply(test1, find_qss, mode = "norm", parms = params.n, 
    iter = 100))
```

```
##    user  system elapsed 
## 4549.44   62.14 4635.06
```

