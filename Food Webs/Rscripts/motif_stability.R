setwd("~/Desktop/GitHub/Ecological Networks/Projects/Food Webs/Rdata")
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
ms.mat <- matrix(unlist(mot.stability), ncol = 13)
colnames(ms.mat) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
ms.mat


save.image("~/Desktop/motstab.Rdata")