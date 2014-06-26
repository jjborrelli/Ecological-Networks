init.path <- "C:/Users/borre_000/Dropbox/Food Web Database/Food_Web/SIZEWEB/" 
setwd(init.path)
sizeweb <- read.csv("sizewebDB.csv", header = T)
dim(sizeweb)
colnames(sizeweb)
length(unique(sizeweb$Web.ID))
