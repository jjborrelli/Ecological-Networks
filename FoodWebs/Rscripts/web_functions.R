###
###
#  Functions for food web structure analyses
###
###

# Function that imports edgelist csv files and outputs graph objects and adjacency matrices
get_webs <- function(directory){  
  # Imports edgelist csv files 
  ## Output is list with two items $graph.list and $adjacency.list
  ### $graph.list is a list of graph objects
  ### $adjacency.list is a list of adjacency matrices
  
  require(igraph)
  require(reshape2)
  
  # Working directory where edgelist csv files are stored
  setwd(directory)
  
  # Get the list of files
  edgelists.files <- as.list(list.files()[grep(".csv",list.files())])
  names(edgelists.files) <- list.files()[grep(".csv",list.files())]
  
  # Get the names of the food webs without the ".csv" extension
  webnames <- unlist(lapply(edgelists.files, strsplit, split = ".", fixed = T))
  webnames <- melt(webnames[seq(1,length(webnames),2)])$value
  # Make all lowercase
  webnames <- tolower(webnames)
  
  # Read in the edgelists
  edgelists <- lapply(edgelists.files, read.csv)
  # Need to switch so that columns are resource - consumer instead of consumer - resource
  for(i in 1:length(edgelists)){
    if(colnames(edgelists[[i]][1]) == "resource"){
      warning("Make sure edgelist csv files have consumer as column 1, and resource as column 2")
    }
    edgelists[[i]]<-edgelists[[i]][,2:1]
  }
  
  # Convert edgelists to graph objects
  edge.graphs <- lapply(edgelists, graph.data.frame)
  names(edge.graphs) <- webnames
  # Convert graphs to adjacency matrices 
  admats <- lapply(edge.graphs, get.adjacency, sparse = FALSE)
  
  return(list(graph.list = edge.graphs, adjacency.list = admats, webnames = webnames))
}

# Function to calculate N, # links, link density, connectance, diameter, average path length, cluster coeff, 
# # loops, and % Top, % Intermediate, % Basal 
get_fw_indices <- function(adj.list, graphs, web){
  require(igraph)
  require(NetIndices)
  require(reshape2)
  
  if(!is.list(adj.list)){
    stop("Input 'adj.list' must be a list of adjacency matrices")
  }
  if(!is.list(graphs)){
    stop("Input 'graphs' should be a list of graph objects")
  }
  if(length(adj.list) != length(graphs)){
    warning("Length of adjacency list and graph list differ")
  }
  # For all food webs calculate general indices
  g.ind <- lapply(adj.list, GenInd)
  # NOTE: Connectance is Links / (N * (N - 1))
  
  # Create a data frame of general indices for all webs
  data <- lapply(g.ind, unlist)
  indices <- data.frame(data[[1]])
  for (i in 2:length(g.ind)){
    indices <- cbind(indices, data[[i]])
  }
  
  # Make indices the columns, webs the rows
  colnames(indices) <- web
  indices <- as.data.frame(t(indices))
  
  
  # Subset the data frame into the desired indices of N, num Links, link density, connectance
  indices <- as.data.frame(indices[,c(1,5,6,7)])
  # Create new column in data frame for web diameter 
  indices <- cbind(indices, D = melt(lapply(graphs,diameter))$value)
  # Add column for average path length
  indices <- cbind(indices, AvPath = melt(lapply(graphs, average.path.length))$value)
  # Add column for "transitivity" or clustering coefficient
  indices <- cbind(indices, Clust = melt(lapply(graphs, transitivity))$value)
  # Loops 
  diagonals <- lapply(adj.list, diag)
  looper2 <- lapply(diagonals, sum)
  # Add column for number of loops in web
  indices$Loop <- unlist(looper2)
  # Calculate proportion Top, Intermediate, Basal
  degrees <- lapply(graphs, degree, mode = "all")
  indegrees <- lapply(graphs, degree, mode = "in")
  outdegrees <- lapply(graphs, degree, mode = "out")
  
  top <- list()
  basal <- list()
  int <- list()
  
  for (i in 1:length(graphs)){
    numBas <- length(indegrees[[i]][which(indegrees[[i]] == 0)])
    numTop <- length(outdegrees[[i]][which(outdegrees[[i]] == 0)])
    basal[[i]] <- (numBas/indices$N[i]) * 100
    top[[i]] <- (numTop/indices$N[i]) * 100
    int[[i]] <- ((indices$N[i] - (numBas + numTop))/indices$N[i]) * 100
  }
  
  indices <- cbind(indices, Top = unlist(top), Int = unlist(int), Bas = unlist(basal))
  
  return(indices)
}

# Function that calculates node properties
get_node_properties <- function(adj.list, web){
  require(igraph)
  require(NetIndices)
  require(reshape2)
  
  if(!is.list(adj.list)){
    stop("Input 'adj.list' must be a list of adjacency matrices")
  }
  
  # For all food webs calculate trophic position and omnivory index
  t.ind <- lapply(adj.list, TrophInd)
  # Convert list of trophic levels to single dataframe
  trophics <- melt(t.ind, id.vars = c("TL", "OI"))
  
  return(trophics)
}

# Count the motif frequencies of a list of graphs
motif_counter <- function(graph.lists, webs = NULL){
  require(igraph)
  
  if(!is.list(graph.lists)){
    stop("The input should be a list of graph objects")
  }
  
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  triad.df <- as.data.frame(triad.matrix)
  
  motif.data.frame <- data.frame(web = webs, s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4, 
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)
  
  return(motif.data.frame)
}

# Create a null distribution of motif frequencies based on rewired webs
null_motifs <- function(graphs, graph.names, sample, iter){
  rewired <- list()
  for(i in 1:sample){
    rewired[[i]] <- lapply(graphs, rewire, mode = "simple", niter = iter)
  }
  rewired.motifs <- lapply(rewired, motif_counter, webs = graph.names)
  nulldis<-do.call(rbind, rewired.motifs)
  return(nulldis)
}

# Generate a niche model food web adjacency matrix
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

###
###
#  Functions for community stability analysis
###
###

# Calculate the eigenvalue with the largest real part of a matrix
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

# Generate a list of niche model webs
niche_maker <- function(n, S, C){
  niche.list <- list()
  for (i in 1:n){
    niche.list[[i]]<- niche.model(S, C)
  }
  return(niche.list)
}

# Iteratively randomize matrix elements, calculate eigenvalue for each iteration
# Returns proportion of eigenvalues with negative real part
test_stability <- function(matrices, rep){
  eigens<-lapply(matrices, function(x){replicate(rep, analyze.eigen(x))})
  qss <- lapply(eigens, function(x){sum(x < 0) / rep})
  return(qss)
}


## Write a function to use on a list of matrices for permutation testing
web_permutation <- function(webmatrices, fixedmar, times, filename, from = 1, to){
  #Input list of adjacency matrices - webmatrices
  #Input whether to fix rows, columns, or both - fixedmar
  #Input number of permuted matrices to generate
  require(vegan)
  if (!is.list(webmatrices)){
    webmatrices <- list(webmatrices)
    warning("not a list")
  }
  len <- length(webmatrices)
  for (i in from:to){
    permuted <- permatfull(webmatrices[[i]], fixedmar = fixedmar, mtype = "prab", times = times)
    permuted.graphs <- lapply(permuted$perm, graph.adjacency)
    p.motif <- motif_counter(permuted.graphs, webs = as.character(1:times))
    file1 <- paste(filename, i, sep = "_")
    write.csv(p.motif, file = paste(file1, "csv", sep = "."))
  }
  #Output is a dataframe of motif frequencies 
}
