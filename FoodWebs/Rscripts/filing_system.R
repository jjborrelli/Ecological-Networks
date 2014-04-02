library(reshape2)
library(plyr)

# create function to break apart the three column format into its component lists
sequencer <- function(x){
 out <- list()
 for (i in 1:nrow(x)){
   if(is.na(x[i, 3])){
     df <- data.frame(x[i,1], x[i,2])
     colnames(df) <- c("consumer", "resource")
     out[[i]] <- df
   }
   if(!is.na(x[i, 3])){
     sequence <- seq(x[i, 2], x[i, 3], 1)
     df <- data.frame(x[i,1],sequence)
     colnames(df) <- c("consumer", "resource")
     out[[i]] <- df
   }
 }
 return(out)
}

# function takes in a three column csv file where the columns are "consumer", "resource" 
# and the third column is the end of a sequence from column 2 to col 3 signifying
# a range of resource nodes, or NA
# output is a two column edgelist "consumer" and "resource"
convert.csv <- function(file){
  data.file <- read.csv(file, header = F)
  colnames(data.file) <- c("consumer", "resource")
  split.data <- split(data.file, data.file$consumer)
  s <- lapply(split.data, sequencer)
  edge <- melt(s, id.vars = c("consumer", "resource"))[,1:2]
  return(edge)
}

#convert all three column csv files from dunne/Peace Lab
setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Food_Web_Edgelists")
edgelists<-as.list(list.files())
readin<-lapply(edgelists, convert.csv)
readin[[1]]
lapply(readin,nrow)
names(readin) <- edgelists

# save dunne files
setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/")

for(i in 1:6){
  write.csv(readin[[i]], file = edgelists[[i]])
}


setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Edgelist")
length(list.files())

# read in iweb matrices and convert to edgelists
setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Food_Web_Edgelists/iweb_matrices/header")
adj.mats <- as.list(list.files())
read.mats <- lapply(adj.mats, read.csv, row.names = 1)
names(read.mats) <- adj.mats
library(igraph)
lapply(read.mats, dim)
make.graphs<- lapply(read.mats, graph.adjacency)
edgelists.iweb <- lapply(make.graphs, get.edgelist)
length(edgelists.iweb)

for (i in 1:length(edgelists.iweb)){
  colnames(edgelists.iweb[[i]]) <- c("resource", "consumer")
}

edgelists.iweb2 <- list()
for (i in 1:length(edgelists.iweb)){
  edgelists.iweb2[[i]] <- matrix(c(edgelists.iweb[[i]][,2], edgelists.iweb[[i]][,1]), ncol = 2)
  colnames(edgelists.iweb2[[i]]) <- c("consumer", "resource")
}


# save the iweb edgelists
setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Edgelist")
for(i in 1:length(edgelists.iweb)){
  write.csv(edgelists.iweb2[[i]], file = adj.mats[[i]], row.names = FALSE)
}

setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Food_Web_Edgelists/iweb_matrices/header")
caribe <- as.list(list.files()[c(8,11,17)])
read.caribe <- lapply(caribe, read.csv, row.names = 1)
names(read.caribe) <- caribe

null.finder <- function(mat){
  empty.rc <- c()
  for (i in 1:dim(mat)[1]){
    if(sum(mat[,i]) == 0){ # Take out columns that sum to zero.
      if(sum(mat[i,]) == 0){empty.rc[i] <- i}
    } 
  }  
  return(empty.rc)
}

found <- lapply(read.caribe, null.finder)
found.vec <- lapply(found, function(x){which(is.na(x) == F)})

zrm.mats <- list()
for(i in 1:3){
  zrm.mats[[i]] <- read.caribe[[i]][-found.vec[[i]], -found.vec[[i]]]
}
zrm.mats<-lapply(zrm.mats,t)
car.g <- lapply(zrm.mats, graph.adjacency)
car.e <- lapply(car.g, get.edgelist)
edgelists.car <- list()
for (i in 1:length(car.e)){
  edgelists.car[[i]] <- matrix(c(car.e[[i]][,2], car.e[[i]][,1]), ncol = 2)
  colnames(edgelists.car[[i]]) <- c("consumer", "resource")
}

# save the iweb edgelists
setwd("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Edgelist")
for(i in 1:3){
  write.csv(edgelists.car[[i]], file = caribe[[i]], row.names = FALSE)
}