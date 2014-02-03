# Bascompte et al. Quantitative network

require(igraph)

setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Quantitative Network/")
bmatrix <- read.csv("interactionsmarine.csv", header = TRUE, row.names = 1)
bmatrix <- as.matrix(bmatrix)
bgraph <- graph.adjacency(bmatrix, mode = "directed", weighted = TRUE)

plot.igraph(bgraph, layout = layout.circle, vertex.size = 1, edge.arrow.size = .25, vertex.label = NA)

bgraph

bel <- get.edgelist(bgraph)
bweights <- E(bgraph)$weight

quantile(bweights)
q75 <- which(bweights > .00239)
q50 <- which(bweights > .00030243)
q25 <- which(bweights > .000033)
q0 <- which(bweights > 0)

b75 <- graph.edgelist(bel[q75,])
b50 <- graph.edgelist(bel[q50,])
b25 <- graph.edgelist(bel[q25,])
ball <- graph.edgelist(bel[q0,])


require(reshape2)
require(devtools)
source_url("https://raw.github.com/jjborrelli/Ecological-Networks/master/Projects/Food%20Webs/Rscripts/web_functions.R")

bsubs <- list(b75, b50, b25, ball)
motif.counts <- motif_counter(bsubs, webs = c("b75", "b50", "b25", "ball"))
null_counts <- null_motifs(bsubs, graph.names = c("b75", "b50", "b25", "ball"), sample = 100, iter = 5000)
null_counts <- split(null_counts, null_counts$web)
#null_counts <- melt(null_counts, id.vars = c("web", "s1", "s2", "s3", "s4", "s5", "d1", "d2", 
#                                             "d3", "d4", "d5", "d6", "d7", "d8"))

nulls <- list(null_counts[[3]][,2:14], null_counts[[2]][,2:14], null_counts[[1]][,2:14], null_counts[[4]][,2:14])
null.mean <- t(sapply(nulls, colMeans))
null.sd <- t(sapply(nulls, FUN = function(x){apply(x, 2, sd)}))
plot(t((motif.counts[,2:14] - null.mean) / null.sd)[,4])

bwebind <- get_fw_indices(adj.list = lapply(bsubs, get.adjacency, sparse = F), 
               graphs = bsubs, web = c("b75", "b50", "b25", "ball"))

bnodeprops <- get_node_properties(adj.list = lapply(bsubs, get.adjacency, sparse = F), 
                                  web = c("b75", "b50", "b25", "ball"))

require(ggplot2)
troplot <- ggplot(bnodeprops, aes(x = TL)) 
troplot <- troplot + geom_histogram(aes(y = ..density..), binwidth = .2, colour = "black", fill = "white") 
troplot + facet_grid(L1 ~ .) + geom_density(alpha = .2, fill = "#FF6666")
