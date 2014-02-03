################################################################################################  
#### Load workspace image                                                                   ####
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Rdata")                   ####
load("chains.Rdata")                                                                        ####
node.props <- read.csv("node properties.csv")                                               ####
################################################################################################
#### Libraries
library(ggplot2)
library(grid)
library(reshape2)
library(igraph)
################################################################################################
#### Plotting trophic level distributions

consumers <- which(round(node.props$TL, 6) >= 2)
tc.plot<-qplot(node.props$TL[consumers], binwidth = .25, geom = "histogram", 
               xlab = "Trophic Position", ylab = "Frequency")
tc.plot <- tc.plot + theme(axis.title.x = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.title.y = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.text.x = element_text(size = 15))
tc.plot <- tc.plot + theme(axis.text.y = element_text(size = 15))
tc.plot + scale_x_continuous(breaks = 1:6)

#setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Figures")
#ggsave(filename = "Trophic_Hist.jpeg", height = 15, width = 15, units = "cm")
#dev.off()

################################################################################################
################################################################################################
#### Plotting chains

graph.chains<-lapply(chain.mats,graph.adjacency)

twospec2<-matrix(c(1,1,
                   2,2),nrow=2,ncol=2,byrow=T)
threespec2<-matrix(c(1,1,
                     3,1,
                     2,2),nrow=3,ncol=2,byrow=T)
fourspec2<-matrix(c(1,1,
                    2,2,
                    0,2,
                    1,3),nrow=4,ncol=2,byrow=T)
fivespec2<-matrix(c(2,1,
                    3,2,
                    1,2,
                    3,3,
                    1,3),nrow=5,ncol=2,byrow=T)
sixspec2<-matrix(c(2,1,
                   3,2,
                   1,2,
                   3,3,
                   1,3,
                   2,4),nrow=6,ncol=2,byrow=T)

layouts<-list(twospec2,threespec2,fourspec2,fivespec2,sixspec2)
text<-c("a","b","c","d","e")

for(i in 1:5){
  E(graph.chains[[i]])$color = "darkslategray4"
  E(graph.chains[[i]], path = c(1:(i+1)))$color = "darkslategrey"
}

#setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Figures")
#jpeg("Chains.jpeg", quality = 100, height = 900, width = 450)

par(mfrow=c(5, 1),mar=c(.5, .5, .5, .5))
for(i in 1:5){
  plot.igraph(graph.chains[[i]], layout = layouts[[i]],
              vertex.size = 40,
              vertex.color = "black",
              vertex.label.color = "white",
              vertex.label.cex = 1,
              edge.width = 3,
              edge.arrow.size = .35,
              frame = T)
  text(2, 1, label = text[i])
}

#dev.off()
################################################################################################
################################################################################################
#### Plotting QSS

qss.df <- as.data.frame(qss.tab)
colnames(qss.df) <- LETTERS[1:7]
qss.df <- melt(qss.df)
qss.df$level <- rep.int(c(1,2,3,4,5,6,7,8,9), 7)

g <- ggplot(qss.df, aes(x = level, y = value))
g <- g + geom_point(aes(shape = qss.df$variable), size = 4)
g <- g + geom_line(aes(linetype = qss.df$variable))
g <- g + theme(legend.title = element_blank()) 
g <- g + theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size = 20))
g <- g + theme(axis.title.x = element_text(size=30))
g <- g + theme(axis.title.y = element_text(size=30))
g <- g + theme(axis.text.x = element_text(size=20))
g <- g + theme(axis.text.y = element_text(size=20))
g <- g + ylab("QSS") + xlab("Chain Length")
g 

#setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Figures")
#ggsave(filename = "QSS_unif.jpeg", height = 15, width = 15, units = "cm")

################################################################################################
################################################################################################
#### Plotting QSS - normal dist

qss.df.n <- as.data.frame(qss.tab.n)
colnames(qss.df.n) <- LETTERS[1:7]
qss.df.n <- melt(qss.df.n)
qss.df.n$level <- rep.int(c(1,2,3,4,5,6,7,8,9), 7)

gn <- ggplot(qss.df.n, aes(x = level, y = value))
gn <- gn + geom_point(aes(shape = qss.df.n$variable), size = 4)
gn <- gn + geom_line(aes(linetype = qss.df.n$variable))
gn <- gn + theme(legend.title = element_blank()) 
gn <- gn + theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size = 20))
gn <- gn + theme(axis.title.x = element_text(size=30))
gn <- gn + theme(axis.title.y = element_text(size=30))
gn <- gn + theme(axis.text.x = element_text(size=20))
gn <- gn + theme(axis.text.y = element_text(size=20))
gn <- gn + ylab("QSS") + xlab("Chain Length")
gn

#setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Figures")
#ggsave(filename = "QSS_norm.jpeg", height = 15, width = 15, units = "cm")
#dev.off()