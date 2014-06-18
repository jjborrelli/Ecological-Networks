Figures for Presentation
========================================================

- Complexity stability

```{r ComplexStability}
require(igraph)

erg1 <- erdos.renyi.game(7, .3, "gnp", directed = T)
erg2 <- erdos.renyi.game(7, .6, "gnp", directed = T)
erg3 <- erdos.renyi.game(7, .9, "gnp", directed = T)

erg4 <- erdos.renyi.game(7, .5, "gnp", directed = T)
erg5 <- erdos.renyi.game(12, .5, "gnp", directed = T)
erg6 <- erdos.renyi.game(17, .5, "gnp", directed = T)

graphs <- list(erg1, erg2, erg3, erg4, erg5, erg6)

par(mfrow = c(2,3), mar = c(.1, .1, .1, .1))
for(i in 1:6){
  plot.igraph(graphs[[i]], frame = T, edge.arrow.size = .25, vertex.label = NA)
  #text(1, 1,tolower(LETTERS[i]))
}

plot.igraph(barabasi.game(30))
```

- Parameter space

  - make a circle:
```{r}
require(ggplot2)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(N1 = xx, N2 = yy))
}

c1 <- circleFun(c(5,5), 5, npoints = 100)
c2 <- circleFun(c(5,5), 1, npoints = 100)

c1 <- cbind(c1, lev = factor(rep("HIGH QSS", 100)))
c2 <- cbind(c2, lev = factor(rep("LOW QSS", 100)))

cir <- rbind(c1, c2)
```

Use brownian motion to simulate stochastically varying populations

```{r}
movement <- function(steps, mode = "normal", p1 = 0, p2 = 1){
  if(mode == "normal"){
    distances <- abs(rnorm(steps, p1, p2))
  }
  if(mode == "lognormal"){
    distances <- rlnorm(steps, p1, p2)
  }
  if(mode == "uniform"){
    distances <- abs(runif(steps, p1, p2))
  }
  return(distances)
}

random_walk.2 <- function(steps, left.p = .5, up.p = .5, mode = "normal", p1 = 0, p2 = 1){
  
  horizontal <- sample(c(-1,1), steps, replace = T, prob = c(left.p, 1-left.p))
  vertical <- sample(c(-1,1), steps, replace = T, prob = c(up.p, 1-up.p))
  moving.horiz <- movement(steps, mode, p1, p2)
  moving.vert <- movement(steps, mode, p1, p2)
  
  walker <- matrix(c(runif(1, 0, 10), runif(1, 0, 10)), nrow = steps+1, ncol = 2, byrow = T)
  for(i in 1:steps){
    # Horizontal Movement
    if(horizontal[i] == 1){
      walker[i+1,1] <- walker[i,1] + moving.horiz[i]
    }
    if(horizontal[i] == -1){
      walker[i+1,1] <- walker[i,1] - moving.horiz[i]
    }
    
    # Vertical Movement
    if(vertical[i] == 1){
      walker[i+1,2] <- walker[i,2] + moving.vert[i]
    }
    if(vertical[i] == -1){
      walker[i+1,2] <- walker[i,2] - moving.vert[i]
    }
  }
  
  colnames(walker) <- c("x", "y")
  walker <- as.data.frame(walker)
  
  return(walker)       
}
        

#w1 <- abs(random_walk.2(25, left.p = .7, up.p = .7))
#w2 <- abs(random_walk.2(25, left.p = .7, up.p = .7))

#w1 <- cbind(w1, lev = factor(rep("HIGH QSS", 26)))
#w2 <- cbind(w2, lev = factor(rep("LOW QSS", 26)))

#walks <- rbind(w1, w2)

```

Create a GIF of different regions of stable parameter spaces to show that higher qss has higher probability of stability. 

```{r}
require(animation)

setwd("C:/Users//borre_000/Documents/Presentations/images/")
saveGIF(
{for(i in c(2:26, 26, 26, 26, 26)){
  g <- ggplot(cir, aes(x = N1, y = N2)) 
  g <- g + geom_polygon(alpha = .5) 
  g <- g + geom_path(data = walks[c(1:i, 27:(i+26)),], aes(x, y)) 
  g <- g + geom_point(data = walks[c(1:i, 27:(i+26)),], aes(x, y), col = "darkgreen", size = 2.5)
  g <- g + scale_x_continuous(limits = c(0, 10)) 
  g <- g + scale_y_continuous(limits = c(0, 10)) 
  g <- g + facet_wrap(~lev)
  print(g)
}}, movie.name = "swissgifEX.gif", interval = .5, nmax = 25, ani.width = 600, ani.height = 600,
outdir = getwd()  
)

```

```{r}
load("C:/Users/borre_000/Desktop/GitHub/Food-Chain-Length/WholeWebSim//chainINFO.Rdata")

a <- which(totalDAT$scenario == c("1/-1", "5/-5", "10/-10"))
b <- which(totalDAT$scenario == c("10/-1", "5/-1", "1/-5", "1/-10"))

ggplot(totalDAT[b,], aes(x = MeanTL, y = QSS)) +
  geom_point(aes(col = ints)) +
  geom_smooth(aes(lty = ints), method = "glm", col = "black") +
  facet_wrap(~scenario)

```

```{r}
m <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:9){
  m[i,i+1] <- 1
}
m[8,10] <- 1
m[5,9] <- 1
m[3,7] <- 1
plot.igraph(graph.adjacency(m), layout = layout.circle)

#---------------------
m <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:9){
  m[i,i+1] <- 1
}
m[1,9] <- 1
m[2,5] <- 1
m[3,c(5,8)] <- 1
m[5,8] <- 1
m[7,10] <- 1
m[8,10] <- 1
m
plot.igraph(graph.adjacency(m), layout = layout.circle)

#---------------------
m <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:9){
  m[i,i+1] <- 1
}
m[1,c(4,9)] <- 1
m[2,c(5,10)] <- 1
m[3,c(5,7,8)] <- 1
m[5,10] <- 1
m[6,8] <- 1
m[7,9] <- 1
m[8,10] <- 1
plot.igraph(graph.adjacency(m), layout = layout.circle)

#---------------------
m <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:9){
  m[i,i+1] <- 1
}
m[1,c(3,7,9)] <- 1
m[2,c(4,5,8,10)] <- 1
m[3,c(6,9)] <- 1
m[4,7] <- 1
m[5,c(7,8,10)] <- 1
m[6,9] <- 1
m[8,10] <- 1
plot.igraph(graph.adjacency(m), layout = layout.circle)

#---------------------
m <- matrix(0, nrow = 10, ncol = 10)
for(i in 1:9){
  m[i,i+1] <- 1
}
m[1,c(3,6,8,10)] <- 1
m[2,c(5,7,8,10)] <- 1
m[3,c(5,7,9)] <- 1
m[4,c(7,8)] <- 1
m[5,c(9,10)] <- 1
m[6,c(8,9)] <- 1
m[7,9] <- 1
m[8,10] <- 1
plot.igraph(graph.adjacency(m), layout = layout.circle)

```