require(bipartite)
require(ggplot2)

# Quantitative pollination networks
setwd("C:/Users/borre_000/Dropbox/Food Web Database/Mutualistic/")
polwebfiles <- list.files("web-of-life_2014-03-12_000716/")[grep("M_PL", list.files("web-of-life_2014-03-12_000716/"))]

polwebs <- list()
for(i in 1:length(polwebfiles)){
  polwebs[[i]] <- read.csv(paste("web-of-life_2014-03-12_000716/",polwebfiles[i], sep = ""), header = T, row.names = 1)
}

pwebs <- lapply(polwebs, as.matrix)
polspecial <- sapply(pwebs, function(x){apply(x, 1, function(y){sum(y>0)/ncol(x)})})
polspecial2 <- sapply(pwebs, function(x){apply(x, 2, function(y){sum(y>0)/nrow(x)})})
lps <- sapply(polspecial, length)
lps2 <- sapply(polspecial2, length)

qpolsp <- sapply(pwebs, function(x){apply(x, 1, function(y){sum(max(y)-y)/(ncol(x)-1)})})
qpolsp2 <- sapply(pwebs, function(x){apply(x, 2, function(y){sum(max(y)-y)/(nrow(x)-1)})})

polrefs <- read.csv("web-of-life_2014-03-12_000716/references.csv")
colnames(polrefs)

df <- data.frame(spec = c(names(unlist(polspecial)), names(unlist(polspecial2))), 
                 special = c(unlist(polspecial), unlist(polspecial2)), 
                 lat = c(rep(polrefs$Latitude, lps),rep(polrefs$Latitude, lps2)), 
                 long = c(rep(polrefs$Longitude, lps),rep(polrefs$Longitude, lps2)), 
                 web = c(as.character(rep(polrefs$ID, lps)),as.character(rep(polrefs$ID, lps2))), 
                 typ = c(as.character(rep(polrefs$Type.of.interactions,lps)),as.character(rep(polrefs$Type.of.interactions,lps2))))

# Quantitative seed disperser networks
sdwebfiles <- list.files("web-of-life_2014-03-12_000816/")[grep("M_SD", list.files("web-of-life_2014-03-12_000816/"))]

sdwebs <- list()
for(i in 1:length(sdwebfiles)){
  sdwebs[[i]] <- read.csv(paste("web-of-life_2014-03-12_000816/",sdwebfiles[i], sep = ""), header = T, row.names = 1)
}

swebs <- lapply(sdwebs, as.matrix)
sspecial <- sapply(swebs, function(x){apply(x, 1, function(y){sum(y>0)/ncol(x)})})
sspecial2 <- sapply(swebs, function(x){apply(x, 2, function(y){sum(y>0)/nrow(x)})})
lss <- sapply(sspecial, length)
lss2 <- sapply(sspecial2, length)

qsdsp <- sapply(swebs, function(x){apply(x, 1, function(y){sum(max(y)-y)/(ncol(x)-1)})})
qsdsp2 <- sapply(swebs, function(x){apply(x, 2, function(y){sum(max(y)-y)/(nrow(x)-1)})})

sdref <- read.csv("web-of-life_2014-03-12_000816/references.csv")

df3 <- data.frame(spec = c(names(unlist(sspecial)),names(unlist(sspecial2))), 
                  special = c(unlist(sspecial),unlist(sspecial2)), 
                  lat = c(rep(sdref$Latitude, lss),rep(sdref$Latitude, lss2)), 
                  long = c(rep(sdref$Longitude, lss),rep(sdref$Longitude, lss2)), 
                  web = c(as.character(rep(sdref$ID, lss)),as.character(rep(sdref$ID, lss2))), 
                  typ = c(as.character(rep(sdref$Type.of.interactions, lss)),as.character(rep(sdref$Type.of.interactions, lss2))))

# Binary plant pollinator networks
polwebfiles.bin <- list.files("web-of-life_2014-03-12_000908/")[grep("M_PL", list.files("web-of-life_2014-03-12_000908/"))]

polwebs.bin <- list()
for(i in 1:length(polwebfiles.bin)){
  polwebs.bin[[i]] <- read.csv(paste("web-of-life_2014-03-12_000908/",polwebfiles.bin[i], sep = ""), header = T, row.names = 1)
}

pwebs.bin <- lapply(polwebs.bin, as.matrix)
polspecial.bin <- sapply(pwebs.bin, function(x){apply(x, 1, function(y){sum(y>0)/ncol(x)})})
polspecial.bin2 <- sapply(pwebs.bin, function(x){apply(x, 2, function(y){sum(y>0)/nrow(x)})})
lpsb <- sapply(polspecial.bin, length)
lpsb2 <- sapply(polspecial.bin2, length)

polref.bin <- read.csv("web-of-life_2014-03-12_000908/references.csv")

df2 <- data.frame(spec = c(names(unlist(polspecial.bin)),names(unlist(polspecial.bin2))),
                  special = c(unlist(polspecial.bin),unlist(polspecial.bin2)), 
                  lat = c(rep(polref.bin$Latitude, lpsb),rep(polref.bin$Latitude, lpsb2)), 
                  long = c(rep(polref.bin$Longitude, lpsb),rep(polref.bin$Longitude, lpsb2)), 
                  web = c(as.character(rep(polref.bin$ID, lpsb)),as.character(rep(polref.bin$ID, lpsb2))), 
                  typ = c(as.character(rep(polref.bin$Type.of.interactions, lpsb)),as.character(rep(polref.bin$Type.of.interactions, lpsb2))))

# Binary seed dispersal networks
swebfiles.bin <- list.files("web-of-life_2014-03-12_000927/")[grep("M_SD", list.files("web-of-life_2014-03-12_000927/"))]

swebs.bin <- list()
for(i in 1:length(swebfiles.bin)){
  swebs.bin[[i]] <- read.csv(paste("web-of-life_2014-03-12_000927/",swebfiles.bin[i], sep = ""), header = T, row.names = 1)
}

swebs.bin <- lapply(swebs.bin, as.matrix)
sspecial.bin <- sapply(swebs.bin, function(x){apply(x, 1, function(y){sum(y>0)/ncol(x)})})
sspecial.bin2 <- sapply(swebs.bin, function(x){apply(x, 2, function(y){sum(y>0)/nrow(x)})})
lssb <- sapply(sspecial.bin, length)
lssb2 <- sapply(sspecial.bin2, length)

sdref.bin <- read.csv("web-of-life_2014-03-12_000927/references.csv")

df4 <- data.frame(spec = c(names(unlist(sspecial.bin)),names(unlist(sspecial.bin2))),
                  special = c(unlist(sspecial.bin),unlist(sspecial.bin2)), 
                  lat = c(rep(sdref.bin$Latitude, lssb),rep(sdref.bin$Latitude, lssb2)),
                  long = c(rep(sdref.bin$Longitude, lssb),rep(sdref.bin$Longitude, lssb2)), 
                  web = c(as.character(rep(sdref.bin$ID, lssb)),as.character(rep(sdref.bin$ID, lssb2))), 
                  typ = c(as.character(rep(sdref.bin$Type.of.interactions, lssb)),as.character(rep(sdref.bin$Type.of.interactions, lssb2))))

sp.dat <- rbind(df, df2, df3, df4)


ggplot(udata, aes(x = abs(lat), y = special)) + geom_point(aes(col = typ)) + geom_smooth(method = "glm", aes(col = typ))

ggplot(sp.dat, aes(x = special, y = ..density..)) + geom_histogram()

aggd <- aggregate(sp.dat$special, list(sp.dat$web), median)
aggd

latitudes <- c(polrefs$Latitude, polref.bin$Latitude, sdref$Latitude, sdref.bin$Latitude)
typ <- c(as.character(polrefs$Type.of.interactions), as.character(polref.bin$Type.of.interactions),
         as.character(sdref$Type.of.interactions), as.character(sdref.bin$Type.of.interactions))
adat <- data.frame(aggd, latitudes, typ)

ggplot(adat, aes(x = abs(latitudes), y = x)) + geom_point(aes(col = typ)) + geom_smooth(method = "glm")

summary(glm(adat$x~adat$latitudes))


udata <- sp.dat[-which(unique(sp.dat$spec) %in% sp.dat$spec),]


qdf <- data.frame(specialization = c(unlist(qpolsp), unlist(qpolsp2), unlist(qsdsp), unlist(qsdsp2)),
                  lat = c(rep(polrefs$Latitude, lps),rep(polrefs$Latitude, lps2), rep(sdref$Latitude, lss),rep(sdref$Latitude, lss2)))

pdi1 <- sapply(pwebs, PDI)
pdi2 <- sapply(swebs, PDI)

qdf <- data.frame(specialization = c(unlist(pdi1), unlist(pdi2)),
                  lat = c(rep(polrefs$Latitude, lps2),rep(sdref$Latitude, lss2)))
ggplot(qdf, aes(x = abs(lat), y = specialization)) + geom_point() + geom_smooth(method = "glm")
