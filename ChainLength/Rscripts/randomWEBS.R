iter = 200
tim <- c()
test <- matrix(nrow = 1000, ncol = iter)
for(k in 1:iter){
  require(NetIndices)
  m <- matrix(0, nrow = 20, ncol = 20)
  m[upper.tri(m)] <- rbinom(45, 1, .4)
  tim[k] <- max(TrophInd(m)$TL)
  for(i in 1:10){
    for(j in 1:10){
      if(m[i,j] == 1){m[j,i] <- -1}
    }
  }
  m
  test[,k] <- eig.analysis(1000, list(m))
  cat(k, "\n")  
}
negs <- apply(test, 2, function(x){sum(x<0)/1000})
negs
plot(tim, negs)

summary(lm(negs~tim))

create_matrices <- function(iter, S, p){
  mLIST <- list()
  diam <- c()
  tim <- c()
  require(NetIndices)
  require(igraph)
  
  for(k in 1:iter){
    if(iter < 1000) cat(k, ", ")
    test <- 0
    while(test < 1){
      m <- matrix(0, nrow = S, ncol = S)
      m[upper.tri(m)] <- rbinom(S*(S-1)/2, 1, p)
      r <- apply(m, 1, function(x){sum(x)>=1})
      c <- apply(m, 2, function(x){sum(x)>=1})
      if(sum(c|r) == S){
        if(is.connected(graph.adjacency(m))){
          tim[k] <- max(TrophInd(m)$TL)
          diam[k] <- diameter(graph.adjacency(m)) 
          for(i in 1:S){
            for(j in 1:S){
              if(m[i,j] == 1){m[j,i] <- -1}
            }
          }
          mLIST[[k]] <- m
          test <- 1 
        }
      } 
    }
  }
  res <- list(mats = mLIST, TL = tim, DI = diam)
  return(res)
}

plot.qss <- function(rando, title){
  require(ggplot2)
  df <- data.frame(TL = rando[[2]], QSS = rando[[3]])
  g <- ggplot(df, aes(x = TL, y = QSS)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste(title))
  print(g)
}

web_qss <- function(iter, S, p, plot = TRUE, title){
  mLIST<- create_matrices(iter, S, p)
  cat("\n", "Matrices Created", "\n", "Getting Eigenvalues")
  test <- eig.analysis(10000, mLIST$mats)
  negs <- apply(test, 2, function(x){sum(x<0)/10000})
  results <- list(mLIST$TL, mLIST$DI, negs)
  if(plot) plot.qss(results, title)
  return(results)
}


p <- seq(.1, .6, by = .025)
results2 <- data.frame(NULL)
system.time(
  for(i in 1:length(p)){
    cat("\n", i, "\n")
    rando <- web_qss(1000, 10, p[i], title = p[i])
    results2 <- rbind(results, data.frame(TL = rando[[1]], EIG = rando[[3]], C = factor(p[i]), DI = rando[[2]]+1))
  }
)

fit <- lm(rando[[2]]~rando[[1]])

require(betareg)
results$EIG[which(results$EIG == 0)] <- 0.000000001
results$EIG[which(results$EIG == 1)] <- 0.999999999
b.mod <- betareg(EIG~TL, data = results)
be <- b.mod$fitted.values


success <- results$EIG * 1000
fail <- 1000 - success 
fit <- glm(matrix(c(success, fail), ncol = 2)~results$TL, family = "binomial")
glm.df <-  data.frame(fit = fit$fitted.values, tl = results$TL)

g <- ggplot(results, aes(x = factor(DI), y = EIG)) + geom_boxplot()#geom_point(aes(col = C)) + ylab("QSS") + xlab("Max Trophic Position")
g <- g + facet_wrap(~C)
g <- g + geom_line(data = data.frame(TL = results$TL, Beta = be), aes(x = TL, y = Beta, col = C), col = "blue")
g <- g + geom_line(data = glm.df, aes(x = tl, y = fit), col = "darkgreen")
g 


plot(fit$fitted.values~results$TL)


ys <- results$EIG[-which(results$EIG <= 0 | results$EIG >= 1)]
nonext <- results$TL[-which(results$EIG <= 0 | results$EIG >= 1)]
beta.fitted <- betareg(ys~nonext, link = "logit")$fitted.values
plot(beta.fitted~nonext)

g + geom_line(data = data.frame(x = nonext, y = beta.fitted), aes(x = x, y = y))


new.res <- results[which(results$C == 0.205263157894737),]

newb <- betareg(EIG~TL, data = new.res)$fitted.values

save.image()
getwd()


m <- matrix(0, nrow = 10, ncol = 10)
m[4:8, 4:8][upper.tri(m[4:8, 4:8])] <- rbinom(10, 1, .3)
m[1:3, 4:10][upper.tri(m[1:3, 4:10])] <- rbinom(15, 1, .3)
m[4:8, 9:10] <- rbinom(10, 1, .3)
m
plot.igraph(graph.adjacency(m))
TrophInd(m)
sum(rowSums(m)<=1)
sum(colSums(m)<=1)
