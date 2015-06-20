library(igraph)
library(NetIndices)
library(scatterplot3d)
library(ggplot2)
library(reshape2)

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


randMTL <- function(x, N, C, niche = T){
  mtl <- c()
  metl <- c()
  diam <- c()
  if(niche){
    for(i in 1:x){
      erg <- niche.model(N, C)
      #ergm <- get.adjacency(erg, sparse = F)
      ti <- TrophInd(erg)$TL
      mtl[i] <- max(ti)
      metl[i] <- mean(ti)
      diam[i] <- diameter(graph.adjacency(erg))
    }
  }else{
    for(i in 1:x){
      erg <- erdos.renyi.game(N, C, "gnp", directed = T)
      ergm <- get.adjacency(erg, sparse = F)
      ti <- TrophInd(ergm)$TL
      mtl[i] <- max(ti)
      metl[i] <- mean(ti)
    }
  }
  return(data.frame(mx = mtl, mu = metl, diam = diam, S = rep(N, x), Con = rep(C, x)))
}

con <- seq(.05, .2, .025)
N = seq(10, 100, 10)

x = 10
rtl <- data.frame(mx = c(), mu = c(), S = c(), Con = c())

system.time(
for(i in 1:length(con)){  
  for(j in 1:length(N)){
    rtl <- rbind(rtl, randMTL(x, N[j], con[i], niche = T)) 
    cat(i, "/", length(con), "-", j, "/", length(N), "\n", sep = "")
  }
}
)

ggplot(rtl, aes(x = Con, y = diam)) + geom_point() + geom_smooth(method = "lm")

dat <- data.frame(S = rtl[,3], Con = rtl[,4], mu = rtl[,3])
meandat <- acast(dat, S~Con, fun.aggregate = mean)
meandat
heatmap(meandat, Rowv = NA, Colv = NA)

dat2 <- aggregate(dat$mu, list(dat$S, dat$Con), mean)
library(car)
Anova(lm(x~Group.1*Group.2, data = dat2), type = "III")
summary(lm(x~Group.1*Group.2, data = dat2))

scatterplot3d(rtl)
scatterplot3d(aggregate(dat$mu, list(dat$S, dat$Con), mean))

