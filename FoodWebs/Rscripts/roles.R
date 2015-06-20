# Design a function to take in a list of graphs and get roles

motif_counter <- function(graph.lists){
  require(igraph)
  
  if(!is.list(graph.lists)){
    stop("The input should be a list of graph objects")
  }
  
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  triad.df <- as.data.frame(triad.matrix)
  
  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4, 
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)
  
  rownames(motif.data.frame) <- names(graph.lists)
  return(motif.data.frame)
}

permutes_rc <- function(mat, iter = 100){
  
  pattern1 <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  pattern2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  count <- 0
  
  mat.list <- list()
  for(i in 1:iter){
    mat.list[[i]] <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  }
  
  while(count < iter){
    srow <- sample(1:nrow(mat), 2)
    scol <- sample(1:ncol(mat), 2)
    
    test <- mat[srow, scol]
    
    if(sum(test == pattern1) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}



listROLES <- function(x){
  lam <- list()
  for(i in 1:nrow(x)){
    lam[[i]] <- x[-i,-i]
  }
  
  qgs <- lapply(lam, graph.adjacency)
  
  m.each <- motif_counter(qgs)
  m.all <- motif_counter(list(graph.adjacency(x)))
  
  m2 <- matrix(nrow = nrow(m.each), ncol = 13)
  for(i in 1:nrow(m.each)){
    m2[i,] <- as.numeric(m.all - m.each[i,])
  }

  return(as.data.frame(m2))
}



getROLES <- function(weblist, bind = F){
  require(igraph)
  require(data.table)
  
  amats <- lapply(weblist, get.adjacency, sparse = F)
  m1 <- lapply(amats, listROLES)
  
  if(bind){
    alldat <- rbindlist(m1)
    return(alldat)
  }else{return(m1)}
}


elist1 <- list(otago.nopar.e, sylt.nopar.e, quick.nopar.e, flens.nopar.e, csm.nopar.e, bsq.nopar.e, epb.nopar.e)
glist1 <- lapply(elist1, graph.data.frame)
alist1 <- lapply(glist1, get.adjacency, sparse = F)

#perm <- lapply(alist1, permutes_rc)

rol <- getROLES(web.graphs)

#perm.g <- lapply(perm, function(x){lapply(x, graph.adjacency)})

#perm.rol <- lapply(perm.g, getROLES)

ltind <- lapply(lapply(web.graphs, get.adjacency, sparse = F), TrophInd)
dftind <- rbindlist(ltind)
dfrol <- rbindlist(rol)
head(data.frame(dfrol, dftind))
heatmap(as.matrix(dfrol))
