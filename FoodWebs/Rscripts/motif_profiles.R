require(cheddar)
require(igraph)
require(NetIndices)

# Get the data for 1984 Tuesday Lake
data(TL84)
# Get the adjacency matrix
tl.mat <- PredationMatrix(TL84)
n1 <- nrow(tl.mat)
# This calculate the motifs after eliminating one species from the network
m84 <- motif_counter(list(graph.adjacency(tl.mat)), web = "TL84")
mot84 <- motif_counter(list(graph.adjacency(tl.mat)), web = "TL84")
for(i in 1:n1){
  temp <- tl.mat[-i, -i]
  name <- paste("TL84", i, sep = "-")
  t.mot <- motif_counter(list(graph.adjacency(temp)), web = name)
  mot84 <- rbind(mot84, t.mot)
}

GenInd(tl.mat)

# Now I want to get the particular set of motifs that each species participates in
sp.profiles <- matrix(nrow = n1, ncol = 13)
colnames(sp.profiles) <- names(mot84[,2:14])
rownames(sp.profiles) <- rownames(tl.mat)
for(i in 1:n1){
  diffs <- mot84[1, 2:14] - mot84[1+i, 2:14]
  sp.profiles[i,] <- as.numeric(diffs)
}


# Get the data for 1986 Tuesday Lake
data(TL86)
# Get the adjacency matrix
tl.mat2 <- PredationMatrix(TL86)
n2 <- nrow(tl.mat2)
# This calculate the motifs after eliminating one species from the network
m86 <- motif_counter(list(graph.adjacency(tl.mat2)), web = "TL86")
mot86 <- motif_counter(list(graph.adjacency(tl.mat2)), web = "TL86")
for(i in 1:n2){
  temp <- tl.mat2[-i, -i]
  name <- paste("TL86", i, sep = "-")
  t.mot <- motif_counter(list(graph.adjacency(temp)), web = name)
  mot86 <- rbind(mot86, t.mot)
}

GenInd(tl.mat2)

# Now I want to get the particular set of motifs that each species participates in
sp.profiles86 <- matrix(nrow = n2, ncol = 13)
colnames(sp.profiles86) <- names(mot86[,2:14])
rownames(sp.profiles86) <- rownames(tl.mat2)
for(i in 1:n2){
  diffs <- mot86[1, 2:14] - mot86[1+i, 2:14]
  sp.profiles86[i,] <- as.numeric(diffs)
}

common <- rownames(tl.mat) %in% rownames(tl.mat2)
common2 <- rownames(tl.mat2) %in% rownames(tl.mat)
common == common2

test <- sp.profiles[sort(rownames(sp.profiles[common,])),] - sp.profiles86[sort(rownames(sp.profiles86[common2,])),]

mot84[,2:14] - mot86[,2:14]
