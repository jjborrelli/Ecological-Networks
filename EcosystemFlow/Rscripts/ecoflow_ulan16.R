setwd("C:/Users/borre_000/Dropbox/Food Web Database/Ecosystem Flow/Ulan")

require(enaR)
require(reshape2)

weblist <- list.files()
names <- unlist(lapply(weblist, strsplit, split = ".", fixed = T))
names <- melt(names[seq(1,length(names),2)])$value
ulan.names <- tolower(names)

networks.list <- list()
for (i in 1:15){
  networks.list[[i]] <- read.scor(weblist[i], type = "edge.list")$flow
}

names(networks.list) <- ulan.names
for(i in 1:15){
  write.csv(networks.list[[i]], file = paste(ulan.names[i], ".csv", sep = ""))
}


ulan.mat <- lapply(networks.list, as.matrix.network)
ulan.struct <- t(sapply(ulan.mat, structure.statistics))
colnames(ulan.struct) <- c("n","L","C","LD","lam1A","mlam1A","lam2A","rho","R","d","no.scc","no.scc.big","pscc")
plot(ulan.struct[,1], ulan.struct[,2])

require(igraph)
ulan.igraph <- lapply(ulan.mat, graph.adjacency)
ulan.motifs <- motif_counter(ulan.mat, webs = ulan.names)
ulan.motifs.null <- null_motifs(ulan.igraph, ulan.names, sample = 500, iter = 200)

websp <- split(ulan.motifs.null, ulan.motifs.null$web)
# Calculate column means and standard deviations for each web
null.web.means <- t(sapply(websp, function(x){
  y<-x[,2:14]
  apply(y, 2, mean)}))
null.web.sd <- t(sapply(websp, function(x){
  y<-x[,2:14]
  apply(y, 2, sd)}))

boxplot((ulan.motifs[,2:14] - null.web.means) / null.web.sd)
abline(h=0)

#######################################
network1 <- read.scor(weblist[1], type = "edge.list")
flows.edge <- network1$flow
hist(flows.edge$value)
quantile(flows.edge$value)
plot.igraph(graph.edgelist(as.matrix(flows.edge[,1:2])), layout = layout.circle, 
            edge.width = (flows.edge$value/max(flows.edge$value)), edge.arrow.size = .25)

flows25q <- flows.edge[which(flows.edge$value>1.39),]
flows50q <- flows.edge[which(flows.edge$value>13.7),]
flows75q <- flows.edge[which(flows.edge$value>63.27),]

plot.igraph(graph.edgelist(as.matrix(flows25q[,1:2])), layout = layout.circle, 
            edge.width = (flows25q$value/10), edge.arrow.size = .25, rescale = F)

plot.igraph(graph.edgelist(as.matrix(flows50q[,1:2])), layout = layout.circle, 
            edge.width = (flows50q$value/10), edge.arrow.size = .25)

plot.igraph(graph.edgelist(as.matrix(flows75q[,1:2])), layout = layout.circle, 
            edge.width = (flows75q$value/10), edge.arrow.size = .25)

graphs <- list(graph.edgelist(as.matrix(flows.edge[,1:2])), graph.edgelist(as.matrix(flows25q[,1:2])),
               graph.edgelist(as.matrix(flows50q[,1:2])), graph.edgelist(as.matrix(flows75q[,1:2])))

test.motifs <- motif_counter(graphs, webs = c("0", "25", "50", "75"))

test.motifs.null <- null_motifs(graphs, graph.names = c("0", "25", "50", "75"), sample = 500, iter = 200)

websp <- split(test.motifs.null, test.motifs.null$web)
# Calculate column means and standard deviations for each web
test.web.means <- t(sapply(websp, function(x){
  y<-x[,2:14]
  apply(y, 2, mean)}))
test.web.sd <- t(sapply(websp, function(x){
  y<-x[,2:14]
  apply(y, 2, sd)}))

boxplot((test.motifs[,2:14] - test.web.means) / test.web.sd)
abline(h=0)



## --------------------------------------------
library(igraph)
library(NetIndices)
library(data.table)
library(ggplot2)

setwd("C:/Users/borre_000/Dropbox/Food Web Database/Ecosystem Flow/Ulan_Edges/")
ulanEDGE <- list()
for(i in 1:length(list.files())){
  ulanEDGE[[i]] <- read.csv(list.files()[i], row.names = 1)
}

q50 <- sapply(ulanEDGE, function(x){quantile(x[,3], .5)})

uEDGE.50 <- list()
for(i in 1:length(ulanEDGE)){
  uEDGE.50[[i]] <- ulanEDGE[[i]][which(ulanEDGE[[i]][,3] >= q50[i]),]
}

uGRAPH <- lapply(uEDGE.50, function(x){graph.edgelist(as.matrix(x[,1:2]))})

uMAT <- lapply(uGRAPH, get.adjacency, sparse = F)

tind <- lapply(uMAT, TrophInd)

n <- sapply(tind, nrow)

webs <- rep(ulan.names, n)

allTIND <- cbind(rbindlist(tind), webs = webs)

hist(allTIND$OI)
qplot(allTIND$TL, allTIND$OI, color = allTIND$webs)

par(mfrow=c(3,5), mar = c(.1, .1, .1, .1))

for(i in 1:15){
  plot(uGRAPH[[i]], layout= layout.fruchterman.reingold, vertex.label = NA, vertex.size = 4, edge.arrow.size = .15)
}

all0 <- sapply(lapply(uGRAPH, degree), function(x){sum(x == 0)})
in0 <- sapply(lapply(uGRAPH, degree, mode = "in"), function(x){sum(x == 0)})
