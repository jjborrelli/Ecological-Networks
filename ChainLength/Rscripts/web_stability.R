require(devtools)
require(igraph)
require(NetIndices)

niche.model<-function(S,C){
  new.mat<-matrix(0,nrow=S,ncol=S)
  ci<-vector()
  niche<-runif(S,0,1)
  r<-rbeta(S,1,((1/(2*C))-1))*niche
  
  for(i in 1:S){
    ci[i]<-runif(1,r[i]/2,niche[i])
  }
  
  r[which(niche==min(niche))]<-.00000001
  
  for(i in 1:S){
    
    for(j in 1:S){
      if(niche[j]>(ci[i]-(.5*r[i])) && niche[j]<(ci[i]+.5*r[i])){
        new.mat[j,i]<-1
      }
    }
  }
  
  new.mat<-new.mat[,order(apply(new.mat,2,sum))]
  return(new.mat)
}


nodes <- seq(10, 100, 10)
connectance <- seq(.05, .45, .025) 

#Create a data frame of node/connectance pairings
   
niche.inputs <- data.frame(N = rep(nodes, each = 17), C = rep(connectance, 10))
 

#Generate 10 niche model food webs for each pairing of node and connectance
   
web_maker <- function(n, x){ # where n is # of webs to make for each pairing and x is a data frame of node/connectance pairings
  # Make sure the column names in the data frame are N and C 
  if(colnames(x[1]) != "N"){
    warning("Check column names")
  }
  if(colnames(x[2]) != "C"){
    warning("Check column names")
  }
  
  web <- list()
  
  for(i in 1:length(x[,1])){
    subweb <- list()
    for(j in 1:n){
      subweb[[j]] <- niche.model(x$N[i], x$C[i])
    }
    web[[i]] <- subweb
  }
  
  return(web)
}
 

#Create a list of niche model food webs for each pair of parameters
   
system.time(
  test1 <- web_maker(10, niche.inputs)
)
 

#Call eigenvalue stability analysis functions from GitHub repository

source_url("https://raw.github.com/jjborrelli/Ecological-Networks/master/Projects/Chain%20Length/Rscripts/chains_functions.R")


#Define the parameters for `find_qss`

params.n <- data.frame(pred1 = 0, pred2 = 2, prey1 = 0, prey2 = 1)


#Analyze the quasi-sign stability of the niche model food webs with `find_qss`

system.time(
  test2 <- lapply(test1, find_qss, mode = "norm", parms = params.n, iter = 100)
)
 