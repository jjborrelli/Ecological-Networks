s.ocean <- read.csv("http://esapubs.org/archive/ecol/E092/097/diet.csv")

colnames(s.ocean)

levels(s.ocean$LOCATION)
m <- split(s.ocean, f = s.ocean$LOCATION)

x <- c()
for (i in 1:228){
  x[i] <- nrow(m[[i]])
}

hist(x, breaks = 25, freq = F, ylim = c(0, .002))
require(MASS)
fitdistr(x, "power")
points(dcauchy(sort(x), 19.59, 21.27), typ = "l", col = "blue")
points(dlnorm(sort(x), 3.26, 1.6996), typ = "l", col = "red")
points(dgamma(sort(x), .43877, .003813), typ = "l", col = "darkslategrey")