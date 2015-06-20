niche.model<-function(S,C){
  require(igraph)
  connected = FALSE
  while(!connected){  
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
    
    connected <- is.connected(graph.adjacency(new.mat))
  }
  return(new.mat)
}


n1 <- niche.model(200, .2)
e1 <- get.edgelist(graph.adjacency(n1))

ro <- sample(1:nrow(e1), 1000)
is.connected(graph.edgelist(cbind(as.character(e1[ro,1]), as.character(e1[ro,2]))))

r <- list()
g <- list()
for(i in 1:10){
  r[[i]] <- sample(1:nrow(e1), 1000)
  g[[i]] <- graph.edgelist(cbind(as.character(e1[r[[i]],1]), as.character(e1[r[[i]],2])))
}
g


source("https://raw.githubusercontent.com/jjborrelli/Ecological-Networks/master/FoodWebs/Rscripts/web_functions.R")

get_fw_indices(adj.list = lapply(g, get.adjacency, sparse = F), graphs = g, web = 1:10)

c1 <- list()
for(i in 1:10){
  c1[[i]] <- which(1:200 %in% sort(as.numeric(colnames(get.adjacency(g[[i]], sparse = F)))) == FALSE)
}
c1 
