library(igraph)
library(NetIndices)
library(reshape2)
library(ggplot2)
library(devtools)
library(vegan)

### Load functions       --------------------------------------------------

url <- "https://raw.github.com/jjborrelli/Ecological-Networks/master/Food%20Webs/Rscripts/web_functions.R"
source_url(url)

### Import food web data  --------------------------------------------------
inputs <- get_webs("~/Dropbox/Food Web Database/Food_Web/Edgelist")
web.graphs <- inputs$graph.list
web.matrices <- inputs$adjacency.list
webnames <- inputs$webnames

### Set Working Directory     --------------------------------------------------
setwd("~/Desktop/GitHub/Ecological-Networks/Food Webs")

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

# Create subset of nodes that are herbivores or higher in trophic position
consumers <- which(round(node.props$TL, 6) >= 2)

# Basic plots of trophic position distributions of subset and whole dataset
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


###  Motifs --------------------------------------------------

motif.df <- motif_counter(edge.graphs, webnames)

#write.csv(motif.df, file = "Tables/motifCOUNTS.csv")
motif.df <- read.table("Tables/motifCOUNTS.csv", header = T, sep = ",", row.names = 1)
sub.counts <- motif.df[,2:14]


### Applying permutation methods -----------------------------------------

p.means <- lapply(permutes, FUN = function(x){colMeans(x[,2:14])})
pmeanmat <- matrix(unlist(p.means), ncol = 13, nrow = 50, byrow = T)
p.sd <- lapply(permutes, FUN = function(x){apply(x[,2:14], 2, sd)})
psdmat <- matrix(unlist(p.sd), ncol = 13, nrow = 50, byrow = T)

zscor <- (mot.mat - pmeanmat) / psdmat
znorm <- zscor/colSums(zscor^2, na.rm = T)
boxplot(znorm)

## Calculate permuted webs and confidence intervals ---------

### Fixed row and column margins ----------------------------
system.time(
  permutes <- web_permutation(web.matrices, fixedmar = "both", times = 1000)
)

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

erg <- list()
for(i in 1:10000){
  erg[[i]] <-  erdos.renyi.game(N[1], L[1], "gnm", directed = TRUE, loops = TRUE)
}

merg <- motif_counter(erg, webs = 1:100)[,2:14]
cmerg <- colMeans(merg)
csdmerg <- apply(merg, 2, sd)
(sub.counts[1,] - cmerg) / csdmerg
sqrt(var(merg[,4])/50)
