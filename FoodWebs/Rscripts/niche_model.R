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

library(igraph)

niche.web<-niche.model(10,.2)

analyze.eigen <- function(m){
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- abs(rnorm(1, 0, 3))
        m[j, i] <- -abs(rnorm(1, 0, 1))
      }
    }
  }
  diag(m) <- -1
  ev <- max(Re(eigen(m)$values))
  return(ev)
}

# generate a list of niche model webs
niche_maker <- function(n, S, C){
  niche.list <- list()
  for (i in 1:n){
    niche.list[[i]]<- niche.model(S, C)
  }
  return(niche.list)
}

test_stability <- function(matrices, rep){
  eigens<-lapply(matrices, function(x){replicate(rep, analyze.eigen(x))})
  qss <- lapply(eigens, function(x){sum(x < 0) / rep})
  return(qss)
}

test.list <- niche_maker(50, 25, .15)
system.time(
test.qss <- test_stability(web.matrices, 1000)
)
qss <- matrix(unlist(test.qss), nrow = 50, ncol = 1)

niche.graphs <- lapply(test.list, graph.adjacency)
niche.motif <- motif_counter(lapply(test.list, graph.adjacency), web = 1:50)
niche.motif.df <- niche.motif[,2:14]
stab.motif <- cbind(qss, motif.df)

fit <- lm(qss ~ s1 + s2 + s3 + s4 + s5 + d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8, data = stab.motif)
summary(fit)
fit1 <- lm(qss ~ s2, data = stab.motif)
summary(fit1)

apply(niche.motif.df, 2, mean)

system.time(
  nullniche <- null_motifs(niche.graphs, graph.names = 1:500, sample = 500, iter = 500)
)

# Split giant dataframe (outputt of null_motifs) into a list with each list item its own web
nichesplit <- split(nullniche, nullniche$web)
# Calculate column means and standard deviations for each web
null.niche.means <- lapply(nichesplit, function(x){
  y<-x[,2:14]
  apply(y, 2, mean)})
null.niche.sd <- lapply(nichesplit, function(x){
  y<-x[,2:14]
  apply(y, 2, sd)})

# Combine all means into a single matrix with each row a web, and each column a motif
# Each value is the average of "sample" number of iterations
null.niche.m <- matrix(unlist(null.niche.means), nrow = 500, ncol = 13, byrow = T)
# Same for standard deviation
null.niche.sd <- matrix(unlist(null.niche.sd), nrow = 500, ncol = 13, byrow = T)
# Redefine the motif frequencies as a matrix
n.mot.mat <- as.matrix(niche.motif.df)

# Calculate the z-score as (freq - mean_rewired) / sd_rewired
n.z.mat <- (n.mot.mat - null.niche.m) / null.niche.sd
n.z.colsum <- colSums(n.z.mat^2)

n.z.mat.norm <- n.z.mat / n.z.colsum
# Boxplot is probably the nicest way to visualize the distributions of z.scores for each motif
boxplot(n.z.mat.norm, ylim = c(-.02,.02))
abline(h = 0)


###############
## TESTING Sim

nw <- niche.model(25, .1)
rownames(nw) <- LETTERS[1:25]
colnames(nw) <- LETTERS[1:25]
tl <- TrophInd(nw)$TL
ornw <- nw[order(tl), order(tl)]
test <- list(nw, ornw)

test_stability(test, 10000)

