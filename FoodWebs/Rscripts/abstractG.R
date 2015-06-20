n <- niche.model(10, .12)
diameter(graph.adjacency(n))
require(igraph)
n1 <- niche_maker(1000, 10, .12)
n2 <- lapply(n1, graph.adjacency)
d <- sapply(n2, diameter)
d
max(d)

gmax <- n2[[which.max(d)]]
gmin <- n2[[which.min(d)]]

library(NetIndices)

tmax <- TrophInd(n1[[which.max(d)]])$TL
tmin <- TrophInd(n1[[which.min(d)]])$TL


E(gmax)$color = "darkslategray"
E(gmax, path = get.diameter(gmax))$color = "blue"
l1 <- matrix(c(sample(1:7,10, replace = T), tmax), ncol = 2)
plot(gmax, layout = l1, edge.arrow.size = .25, margin = 0)


E(gmin)$color = "darkslategray"
E(gmin, path = get.diameter(gmin))$color = "blue"
l2 <- matrix(c(1,2,3,4,5,1,2,3,4,5, tmin), ncol = 2)
plot(gmin, layout = l2, edge.arrow.size = .25, margin = 0)
