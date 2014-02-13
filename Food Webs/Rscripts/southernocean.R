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


tent <- data.frame(pred = m[[211]]$PREDATOR_COMMON_NAME, prey = m[[211]]$PREY_COMMON_NAME,
                 frac = m[[211]]$FRACTION_OCCURRENCE)

tent.g <- graph.edgelist(as.matrix(tent[,1:2]))


plot.igraph(tent.g, layout = layout.circle)

m[[4]]
s.sand <- data.frame(pred = m[[4]]$PREDATOR_NAME, prey = m[[4]]$PREY_NAME,
                   frac = m[[4]]$FRACTION_OCCURRENCE)

s.sand.g <- graph.edgelist(unique(as.matrix(s.sand[,1:2])))


plot.igraph(s.sand.g, layout = layout.circle)


# List of graphs by location -----------------------
m <- split(s.ocean, f = s.ocean$LOCATION)  
location.g <- list()  
for (i in 1:length(levels(s.ocean$LOCATION))){
  
  el.df <- data.frame(pred = m[[i]]$PREDATOR_NAME, prey = m[[i]]$PREY_NAME)
  
  g <- graph.edgelist(unique(as.matrix(el.df[,1:2])))
  
  location.g[[i]] <- g
  
  print(i)
}

location.g[[4]]

so.ode <- as.character(s.ocean$OBSERVATION_DATE_END)
so.ode.split <- strsplit(so.ode, split = "/")
length(so.ode.split)

year <- c()
for(i in 1:length(so.ode.split)){
  year[i] <- so.ode.split[[i]][3]
}
head(year)
s.ocean2 <- cbind(s.ocean, year)

m2 <- split(s.ocean2, f = s.ocean2$year)
levels(s.ocean2$year)

x2 <- c()
for (i in 1:length(levels(s.ocean2$year))){
  x2[i] <- nrow(m2[[i]])
}


year <- c()
for(i in 1:length(so.ode.split)){
  year[i] <- so.ode.split[[i]][3]
}
s.ocean2 <- cbind(s.ocean, year)

m2 <- split(s.ocean2, f = s.ocean2$year)
year.g <- list()  
for (i in 1:length(levels(s.ocean2$year))){
  
  el.df <- data.frame(pred = m[[i]]$PREDATOR_NAME, prey = m[[i]]$PREY_NAME)
  
  g <- graph.edgelist(as.matrix(el.df[,1:2]))
  
  year.g[[i]] <- g
  
  print(i)
}

year.g
