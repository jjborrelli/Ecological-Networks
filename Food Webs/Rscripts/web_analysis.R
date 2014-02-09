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
fw.indices <- read.csv("Tables/FWindices.csv")

# Calculate node properties of each food web  ------------------------------
#node.props <- get_node_properties(web.matrices, webnames)

# Save node.props as a csv files
#write.csv(node.props, file = "Tables/NODEproperties.csv")
node.props <- read.csv("Tables/NODEproperties.csv")

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
motif.df <- read.table("Tables/motifCOUNTS.csv", header = T, sep = ",")

motif.means <- apply(motif.df[,3:15], 2, mean)
motif.sd <- apply(motif.df, 2, sd)

motif.df.mc <- motif.df - motif.means         # Mean centered dataframe

motif.quantiles <- apply(z.mat.norm, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

# This function generates a list of lists of each of the graphs with rewired edges (preserving degree distribution)
# It then loops through the list with the motif counter function created above
# The resulting list of dataframes (n = sample) is compiled into one large dataframe of length 50n  

# Here I am running the function on my 50 web list (web.graphs) with names = webnames
# sample is the number of times I want to run the rewire function, iter is the number of rewiring iterations
system.time(
nulltest <- null_motifs(web.graphs, graph.names = webnames, sample = 1000, iter = 1000)
)
## This took 105 minutes, rep = 1500 
## This took 22 minutes, rep = 1000

nullmeans <- apply(nulltest[,2:14], 2, mean)
nullsd <- apply(nulltest[,2:14], 2, sd)
plot((motif.means - nullmeans) / nullsd, typ = "o", lty = 2, ylim = c(-.5, .5), xlim = c(1, 13))
abline(h = 0)

# Split giant dataframe (outputt of null_motifs) into a list with each list item its own web
websplit <- split(nulltest, nulltest$web)
# Calculate column means and standard deviations for each web
null.web.means <- lapply(websplit, function(x){
  y<-x[,2:14]
  apply(y, 2, mean)})
null.web.sd <- lapply(websplit, function(x){
  y<-x[,2:14]
  apply(y, 2, sd)})

# Combine all means into a single matrix with each row a web, and each column a motif
# Each value is the average of "sample" number of iterations
nwm <- matrix(unlist(null.web.means), nrow = 50, ncol = 13, byrow = T)
# Same for standard deviation
nwsd <- matrix(unlist(null.web.sd), nrow = 50, ncol = 13, byrow = T)
# Redefine the motif frequencies as a matrix
mot.mat <- as.matrix(motif.df[,3:15])

# Calculate the z-score as (freq - mean_rewired) / sd_rewired
z.mat <- (mot.mat - nwm) / nwsd

z.mat.norm <- z.mat / colSums(z.mat^2)

# Boxplot is probably the nicest way to visualize the distributions of z.scores for each motif
par(mfrow = c(1,2))

zmat.df <- melt(as.data.frame(z.mat))
mot.pro <- ggplot(zmat.df, aes(x = variable, y = value)) + geom_boxplot()
mot.pro + geom_hline(aes(yintercept = 0))
ggsave(filename = "motif_profile1.jpeg", height = 15, width = 15, units = "cm")

zmatnorm.df <- melt(as.data.frame(z.mat.norm))
mot.pro.norm <- ggplot(zmatnorm.df, aes(x = variable, y = value), na.rm = T) + geom_boxplot()
mot.pro.norm <- mot.pro.norm + geom_hline(aes(yintercept = 0)) 
mot.pro.norm <- mot.pro.norm + geom_point(aes(x = 1:13, y = motif.quantiles[1,]), color = "blue")
mot.pro.norm + geom_point(aes(x = 1:13, y = motif.quantiles[2,]), color = "blue")
#ggsave(filename = "motif_profile2.jpeg", height = 15, width = 15, units = "cm")

###########
## Allesina and Pascual 2008 data
apwebs <- c("benguela", "bridgebrook", "broom", "chesapeake", "littlerock", "reef", "shelf", "stmarks", "stmartin")
ap.motifs <- z.mat[c(3, 6, 8, 13, 31, 38, 39, 41,42),]
rownames(ap.motifs) <- apwebs
ap.qss <- matrix(c(.79, .86, .96, .97, .80, .75, .72, .89, .87), ncol = 1)
ap.motif.qss <- as.data.frame(cbind(ap.motifs, ap.qss))

fittest <- lm(ap.qss ~ s1 + s2 + s3 + s4 + s5, data = ap.motif.qss)
summary(fittest)

nulldistr.web <- do.call(rbind, websplit)


###  Trying different null models ---------------------------------

# Playing with the permatfull function
?permatfull
testing <- web.matrices[[1]]

# Test the permatfull function ----------------------------

permuted2 <- permatfull(testing, fixedmar = "both", mtype = "prab", times = 1000)
plot(permuted)

permuted.graphs2 <- lapply(permuted2$perm, graph.adjacency)

p.motif2 <- motif_counter(permuted.graphs2, webs = as.character(1:1000))
pmeans <- colMeans(p.motif2[2:14])
boxplot(p.motif2[2:14])
quant <- apply(p.motif2[2:14], 2, quantile, probs = c(0.975, 0.025))

plot(pmeans, ylim = c(-1, 1500))
points(quant[1,], col = "blue")
points(quant[2,], col = "blue")
points(tes, col = "red")

### Applying permutation methods -----------------------------------------

p.means <- lapply(permutes, FUN = function(x){colMeans(x[,2:14])})
pmeanmat <- matrix(unlist(p.means), ncol = 13, nrow = 50, byrow = T)
p.sd <- lapply(permutes, FUN = function(x){apply(x[,2:14], 2, sd)})
psdmat <- matrix(unlist(p.sd), ncol = 13, nrow = 50, byrow = T)

zscor <- (mot.mat - pmeanmat) / psdmat
znorm <- zscor/colSums(zscor^2, na.rm = T)
boxplot(znorm)

## Calculate permuted webs and confidence intervals ---------
system.time(
  permutes <- web_permutation(web.matrices, fixedmar = "both", times = 1000)
)

permint.both<- sapply(permutes, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.both, file = "Tables/permutedCI_both.csv")

system.time(
  permutes.row <- web_permutation(web.matrices, fixedmar = "rows", times = 1000)
)

permint.row <- sapply(permutes.row, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.row, file = "Tables/permutedCI_row.csv")

system.time(
  permutes.col <- web_permutation(web.matrices, fixedmar = "columns", times = 1000)
)

permint.col <- sapply(permutes.col, FUN = function(x){apply(x[,2:14], 2, quantile, probs = c(0.975, 0.025))})
#write.csv(permint.col, file = "Tables/permutedCI_col.csv")