install.packages("networkD3")
library(networkD3)
# Create fake data
src <- c("A", "A", "A", "A",
         "B", "B", "C", "C", "D")
target <- c("B", "C", "D", "J",
            "E", "F", "G", "H", "I")
networkData <- data.frame(src, target)

# Plot
simpleNetwork(networkData)


?forceNetwork
library(igraph)
forceNetwork(data.frame(X1=LETTERS[id14[,1]], X2=LETTERS[id14[,2]]), Source = "X1", Target = "X2", Nodes = data.frame(N = LETTERS[1:4], G = c(1,1,1,1)), NodeID = "N", Group = "G")

simpleNetwork(data.frame(id14))



# hive plots
install.packages("HiveR")
library(HiveR)
library(igraph)

erg <- erdos.renyi.game(100, .2, "gnp")
elerg <- get.edgelist(erg)

H <- edge2HPD(data.frame(elerg))
plotHive(H)
