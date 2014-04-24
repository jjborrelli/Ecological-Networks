library(igraph)
library(NetIndices)
library(reshape2)
library(ggplot2)
library(devtools)
library(vegan)

### Load functions       --------------------------------------------------

url <- "https://raw.githubusercontent.com/jjborrelli/Ecological-Networks/master/FoodWebs/Rscripts/web_functions.R"
source_url(url)

source("~/Desktop/GitHub/Ecological-Networks/FoodWebs/Rscripts/web_functions.R")

### Import food web data  --------------------------------------------------
inputs <- get_webs("~/Dropbox/Food Web Database/Food_Web/Edgelist")
web.graphs <- inputs$graph.list
web.matrices <- inputs$adjacency.list
webnames <- inputs$webnames

### Set Working Directory     --------------------------------------------------
setwd("~/Desktop/GitHub/Ecological-Networks/FoodWebs")

### Calculate common food web indices/statistics  ----------------------------
#fw.indices <- get_fw_indices(web.matrices, web.graphs, webnames)
#fw.indices[1:2,]
#write.csv(fw.indices, file = "Tables/FWindices.csv")
fw.indices <- read.csv("Tables/FWindices.csv", row.names = 1)

# Calculate node properties of each food web  ------------------------------
#node.props <- get_node_properties(web.matrices, webnames)

# Save node.props as a csv files
#write.csv(node.props, file = "Tables/NODEproperties.csv")
node.props <- read.csv("Tables/NODEproperties.csv", row.names = 1)
between <- melt(lapply(web.graphs, betweenness))

maxTL <- aggregate(node.props$TL, by = list(node.props$L1), max)

node.props <- cbind(betweenness = between[,1], node.props)
# Create subset of nodes that are herbivores or higher in trophic position
consumers <- which(round(node.props$TL, 6) >= 2)

# Basic plots of trophic position distributions of subset and whole dataset
hist(aggregate(node.props$TL, list(node.props$L1), max)[,2], breaks = seq(2,7,.25))

hist(node.props$TL[consumers], freq = F, breaks = 20, xlim = c(2, 6))

hist(node.props$TL, breaks = 40, freq =F)

plot(node.props$OI ~ node.props$TL)

# ggplot of distribution of trophic positions equal or higher than 2
tc.plot<-qplot(node.props$TL[consumers], binwidth = .8, geom = "histogram", 
             xlab = "Trophic Position", ylab = "Frequency")
tc.plot <- tc.plot + theme(axis.title.x = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.title.y = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.text.x = element_text(size = 15))
tc.plot <- tc.plot + theme(axis.text.y = element_text(size = 15))
tc.plot #+ scale_x_continuous(breaks = 1:6)

# For saving the image to desktop, or working directory
#setwd("~/Desktop") 
#ggsave(filename = "Trophic_Distribution_3.jpeg", height = 15, width = 15, units = "cm")
#dev.off()

### Some plots to explore fw stats
plot(Ltot ~ N, data = fw.indices)
abline(lm(Ltot ~ N, data = fw.indices))
summary(lm(Ltot ~ N, data = fw.indices))

diam.plot<-qplot(fw.indices$D, binwidth = .5, geom = "histogram", 
               xlab = "Diameter", ylab = "Frequency")
diam.plot <- diam.plot + theme(axis.title.x = element_text(size=20))
diam.plot <- diam.plot + theme(axis.title.y = element_blank())
diam.plot <- diam.plot + theme(axis.text.x = element_text(size=15))
diam.plot <- diam.plot + theme(axis.text.y = element_text(size=15))
diam.plot + scale_x_continuous(breaks = 0:9)

#setwd("~/Desktop")
#ggsave(filename = "DaimDistribution2.jpeg", height = 15, width = 15, units = "cm")
#dev.off()


#### Mean center and standardize the data frame

fw.indices
fw.indices.mean <- apply(fw.indices, 2, mean)
fw.indices.sd <- apply(fw.indices, 2, sd)

fw.indices.mc <- fw.indices - fw.indices.mean
fw.indices.stand <- fw.indices.mc / fw.indices.sd

pca.fw <- princomp(fw.indices)
loadings(pca.fw)
plot(predict(pca.fw)[,1:2])
biplot(pca.fw)


###  Motifs --------------------------------------------------

#motif.df <- motif_counter(edge.graphs, webnames)

#write.csv(motif.df, file = "Tables/motifCOUNTS.csv")
motif.df <- read.table("Tables/motifCOUNTS.csv", header = T, sep = ",", row.names = 1)
sub.counts <- motif.df[,2:14]
row.names(sub.counts) <- motif.df[,1]
alt.sub.counts <- sub.counts[-c(12, 18, 28),]
pca.mot <- princomp(sub.counts)
loadings(pca.mot)
plot(predict(pca.mot)[,1:2])

biplot(pca.mot, cex = c(.5, .75), pc.biplot = T)

### Applying permutation methods -----------------------------------------

p.means <- lapply(permutes, FUN = function(x){colMeans(x[,2:14])})
pmeanmat <- matrix(unlist(p.means), ncol = 13, nrow = 50, byrow = T)
p.sd <- lapply(permutes, FUN = function(x){apply(x[,2:14], 2, sd)})
psdmat <- matrix(unlist(p.sd), ncol = 13, nrow = 50, byrow = T)

zscor <- (mot.mat - pmeanmat) / psdmat
znorm <- zscor/colSums(zscor^2, na.rm = T)
boxplot(znorm)

## Calculate permuted webs and confidence intervals ---------

### Fixed row and column sums ----------------------------
system.time(
  permutes <- web_permutation(web.matrices, fixedmar = "both", times = 1000)
)

permean.both<- sapply(permutes, FUN = function(x){apply(x[,2:14], 2, mean)})
persd.both<- sapply(permutes, FUN = function(x){apply(x[,2:14], 2, sd)})

z.both <- (sub.counts - t(permean.both)) / t(persd.both)
z.both[is.nan(as.matrix(z.both))] <- 0
z.norm <- apply(z.stand, 2, FUN = function(x){(x * abs(sum(x)))/sqrt(sum(x^2))})
#write.csv(z.norm, file = "Tables/zscore_both.csv")
z.norm <- read.csv("Tables/zscore_both.csv", row.names = 1)
zeros <- which(as.numeric(rowSums(sub.counts[,6:13])) == 0)
boxplot(z.norm)
abline(h = 0)



permint.both<- sapply(permutes, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.both, file = "Tables/permutedCI_both.csv")
perm.both <- read.csv("Tables/permutedCI_both.csv", row.names = 1)

upper.quant.both <- t(perm.both[odds,])
lower.quant.both <- t(perm.both[evens,])

colSums(motif.df[,2:14] >= upper.quant.both) / 50
colSums(motif.df[,2:14] <= lower.quant.both) / 50

## 
## Calculate mean and sd for permuted matrices
pb.means <- t(sapply(permutes, FUN = function(x){apply(x[,2:14], 2, mean)}))
pb.sd <- t(sapply(permutes, FUN = function(x){apply(x[,2:14], 2, sd)}))
##
## Calculate z scores
z.b <- (sub.counts - pb.means) / pb.sd 

##
## Calculate standard error
se.zb <- apply(z.b, 2, FUN = function(x){sqrt(var(x[!is.na(x)])/length(x))})
z.bmean <- colMeans(z.b, na.rm = T)

z.high.b <- z.bmean + se.zb
z.low.b <- z.bmean - se.zb

##
## Plot effect sizes 
plot(z.bmean, ylim = c(-5, 15))
points(z.high.b, col = "blue", pch = 16)
points(z.low.b, col = "blue", pch = 16)
abline(h = 0, lty = 2)

### Fixed row margins ---------------------------------------
system.time(
  permutes.row <- web_permutation(web.matrices, fixedmar = "rows", times = 1000)
)

permint.row <- sapply(permutes.row, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.row, file = "Tables/permutedCI_row.csv")
perm.row <- read.csv("Tables/permutedCI_row.csv", row.names = 1)

upper.quant.row <- t(perm.row[odds,])
lower.quant.row <- t(perm.row[evens,])

colSums(motif.df[,2:14] >= upper.quant.row) / 50
colSums(motif.df[,2:14] <= lower.quant.row) / 50

## 
## Calculate mean and sd for permuted matrices
pr.means <- t(sapply(permutes.row, FUN = function(x){apply(x[,2:14], 2, mean)}))
pr.sd <- t(sapply(permutes.row, FUN = function(x){apply(x[,2:14], 2, sd)}))
##
## Calculate z scores
z.r <- (sub.counts - pr.means) / pr.sd 
z.r[is.nan(as.matrix(z.r))] <- 0
zr.stand <- apply(z.r, 2, FUN = function(x){x/abs(sum(x))})
#write.csv(zr.stand, file = "Tables/zscore_row.csv")
boxplot(zr.stand)
abline(h = 0)
apply(zr.stand,2, FUN = function(x){sum(x>0)/length(x)})


##
## Calculate standard error
se.zr <- apply(z.r, 2, FUN = function(x){sqrt(var(x[!is.na(x)])/length(x))})
z.rmean <- colMeans(z.r, na.rm = T)

z.high.r <- z.rmean + se.zr
z.low.r <- z.rmean - se.zr

##
## Plot effect sizes 
plot(z.rmean, ylim = c(-10, 10))
points(z.high.r, col = "blue", pch = 16)
points(z.low.r, col = "blue", pch = 16)
abline(h = 0, lty = 2)

### Fixed column margins ------------------------------------
system.time(
  permutes.col <- web_permutation(web.matrices, fixedmar = "columns", times = 1000)
)


permint.col <- sapply(permutes.col, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.col, file = "Tables/permutedCI_col.csv")
perm.col <- read.csv("Tables/permutedCI_col.csv", row.names = 1)

upper.quant.col <- t(perm.col[odds,])
lower.quant.col <- t(perm.col[evens,])

colSums(motif.df[,2:14] >= upper.quant.col) / 50
colSums(motif.df[,2:14] <= lower.quant.col) / 50

## 
## Calculate mean and sd for permuted matrices
pc.means <- t(sapply(permutes.col, FUN = function(x){apply(x[,2:14], 2, mean)}))
pc.sd <- t(sapply(permutes.col, FUN = function(x){apply(x[,2:14], 2, sd)}))
##
## Calculate z scores
z.c <- (sub.counts - pc.means) / pc.sd 
z.c[is.nan(as.matrix(z.c))] <- 0
zc.stand <- apply(z.c, 2, FUN = function(x){x/abs(sum(x))})
#write.csv(zc.stand, file = "Tables/zscore_col.csv")
boxplot(zc.stand)
abline(h = 0)
apply(zc.stand,2, FUN = function(x){sum(x>0)/length(x)})


##
## Calculate standard error
se.zc <- apply(z.c, 2, FUN = function(x){sqrt(var(x[!is.na(x)])/length(x))})
z.cmean <- colMeans(z.c, na.rm = T)

z.high.c <- z.cmean + se.zc
z.low.c <- z.cmean - se.zc

##
## Plot effect sizes 
plot(z.cmean, ylim = c(-10, 10))
points(z.high.c, col = "blue", pch = 16)
points(z.low.c, col = "blue", pch = 16)
abline(h = 0, lty = 2)

### Applying null model with random webs ------------------------------

N <- fw.indices$N
L <- fw.indices$Ltot

list_erg <- function(n, l, times = 1000, loops = TRUE){
  erg <- list()
  for(i in 1:times){
    erg[[i]] <-  erdos.renyi.game(n, l, "gnm", directed = TRUE, loops = loops)
  } 
  return(erg)
}

erg_counter <- function(N, L, iter = 1000, loop = TRUE, webs = "NONE"){
   if(!length(N) == length(L)){
    stop("N and L vectors not of same length")
   }
   if(webs == "NONE"){
    webs <- 1:iter
  }
  
  mot.ergs <- list()
  for(fw in 1:length(N)){
    ergs <- list_erg(N[fw], L[fw], times = iter, loops = loop)
    mot.ergs[[fw]] <- motif_counter(ergs, webs = webs)[,2:14]
  }
  return(mot.ergs)
}

testerg <- erg_counter(N, L, iter = 100, loop = TRUE)


