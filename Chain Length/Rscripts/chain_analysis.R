################################################################################################
################################################################################################
#### Load source functions                                                                  ####
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Rscripts")                ####
source("chains_functions.R")                                                                ####
################################################################################################  
#### Load workspace image                                                                   ####
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Rdata")                   ####
load("chains.Rdata")                                                                        ####
################################################################################################
################################################################################################

chain.mats <- list_chains(2,10)

##### With Unform distribution
params.u <- data.frame(pred1 = c(0, 0, 0, 0, 0, 0, 0), pred2 = c(10, 10, 10, 10, 5, 3, 1), 
                     prey1 = c(-5, -1, -.1, -.01, -1, -1, -1), prey2 = c(0, 0, 0, 0, 0, 0, 0))

qss.tab <- find_qss(chain.mats, mode = "unif", parms = params.u, iter = 10000)

##### With Normal Distribution
params.n <- data.frame(pred1 = c(0, 0, 0, 0, 0, 0, 0), pred2 = c(1, 2, 3, 4, 5, 6, 7), 
                       prey1 = c(0, 0, 0, 0, 0, 0, 0), prey2 = c(1, 1, 1, 1, 1, 1, 1))

qss.tab.n <- find_qss(chain.mats, mode = "norm", parms = params.n, iter = 10000)

################################################################################################
################################################################################################
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Rdata")
save.image("chains.Rdata")