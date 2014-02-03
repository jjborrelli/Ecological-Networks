library(igraph)
library(NetIndices)
library(reshape2)
library(ggplot2)
library(devtools)


### Restore workspace    -------------------------------------------------
                                                                         
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Food Webs/Rdata")   
load("foodweb.Rdata")                                                    

### Load functions       --------------------------------------------------

url <- "https://raw.github.com/jjborrelli/Ecological-Networks/master/Projects/Food%20Webs/Rscripts/web_functions.R"
source_url(url)

### Import food web data  --------------------------------------------------
inputs <- get_webs("~/Desktop/GitHub/Ecological Networks/Database/Food_Web/Edgelist")
web.graphs <- inputs$graph.list
web.matrices <- inputs$adjacency.list
webnames <- inputs$webnames

### Calculate common food web indices/statistics  ----------------------------
fw.indices <- get_fw_indices(web.matrices, web.graphs, webnames)
fw.indices[1:2,]

# Calculate node properties of each food web  ------------------------------
node.props <- get_node_properties(web.matrices, webnames)

# Save node.props as a csv files
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Chain Length/Rdata")
write.csv(node.props, "node properties.csv")

# Create subset of nodes that are herbivores or higher in trophic position
consumers <- which(round(node.props$TL, 6) >= 2)

# Basic plots of trophic position distributions of subset and whole dataset
hist(node.props$TL[consumers], freq = F, breaks = 20, xlim = c(2, 6))

hist(node.props$TL, breaks = 40, freq =F)

plot(node.props$OI ~ node.props$TL)

# ggplot of distribution of trophic positions equal or higher than 2
tc.plot<-qplot(node.props$TL[consumers], binwidth = .25, geom = "histogram", 
             xlab = "Trophic Position", ylab = "Frequency")
tc.plot <- tc.plot + theme(axis.title.x = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.title.y = element_text(size = 20))
tc.plot <- tc.plot + theme(axis.text.x = element_text(size = 15))
tc.plot <- tc.plot + theme(axis.text.y = element_text(size = 15))
tc.plot + scale_x_continuous(breaks = 1:6)

# For saving the image to desktop, or working directory
#setwd("~/Desktop") 
#ggsave(filename = "Trophic_Distribution_3.jpeg", height = 15, width = 15, units = "cm")
#dev.off()

### Some plots to explore fw stats
plot(Ltot ~ N, data = fw.indices)
abline(lm(Ltot ~ N, data = fw.indices))
summary(lm(Ltot ~ N, data = fw.indices))

diam.plot<-qplot(fw.indices$D, binwidth = .5, geom = "histogram", 
               xlab = "Diameter", ylab = "Frequency")
diam.plot <- diam.plot + theme(axis.title.x = element_text(size=20))
diam.plot <- diam.plot + theme(axis.title.y = element_blank())
diam.plot <- diam.plot + theme(axis.text.x = element_text(size=15))
diam.plot <- diam.plot + theme(axis.text.y = element_text(size=15))
diam.plot + scale_x_continuous(breaks = 0:9)

#setwd("~/Desktop")
#ggsave(filename = "DaimDistribution2.jpeg", height = 15, width = 15, units = "cm")
#dev.off()


#### Mean center and standardize the data frame

fw.indices
fw.indices.mean <- apply(fw.indices, 2, mean)
fw.indices.sd <- apply(fw.indices, 2, sd)

fw.indices.mc <- fw.indices - fw.indices.mean
fw.indices.stand <- fw.indices.mc / fw.indices.sd

pca.fw <- princomp(fw.indices)
loadings(pca.fw)
plot(predict(pca.fw)[,1:2])


###  Motifs --------------------------------------------------

motif.df <- motif_counter(edge.graphs, webnames)

motif.means <- apply(motif.df, 2, mean)
motif.sd <- apply(motif.df, 2, sd)

motif.df.mc <- motif.df - motif.means         # Mean centered dataframe

# This function generates a list of lists of each of the graphs with rewired edges (preserving degree distribution)
# It then loops through the list with the motif counter function created above
# The resulting list of dataframes (n = sample) is compiled into one large dataframe of length 50n  

# Here I am running the function on my 50 web list (web.graphs) with names = webnames
# sample is the number of times I want to run the rewire function, iter is the number of rewiring iterations
system.time(
nulltest <- null_motifs(web.graphs, graph.names = webnames, sample = 1500, iter = 1000)
)
## This took 105 minutes

nullmeans <- apply(nulltest[,2:14], 2, mean)
nullsd <- apply(nulltest[,2:14], 2, sd)
plot((motif.means - nullmeans) / nullsd, typ = "o", lty = 2, ylim = c(-.5, .5), xlim = c(1, 13))
abline(h = 0)

# Split giant dataframe (outputt of null_motifs) into a list with each list item its own web
websplit <- split(nulltest, nulltest$web)
# Calculate column means and standard deviations for each web
null.web.means <- lapply(websplit, function(x){
  y<-x[,2:14]
  apply(y, 2, mean)})
null.web.sd <- lapply(websplit, function(x){
  y<-x[,2:14]
  apply(y, 2, sd)})

# Combine all means into a single matrix with each row a web, and each column a motif
# Each value is the average of "sample" number of iterations
nwm <- matrix(unlist(null.web.means), nrow = 50, ncol = 13, byrow = T)
# Same for standard deviation
nwsd <- matrix(unlist(null.web.sd), nrow = 50, ncol = 13, byrow = T)
# Redefine the motif frequencies as a matrix
mot.mat <- as.matrix(motif.df[,2:14])

# Calculate the z-score as (freq - mean_rewired) / sd_rewired
z.mat <- (mot.mat - nwm) / nwsd

z.mat.norm <- z.mat / colSums(z.mat^2)

# Boxplot is probably the nicest way to visualize the distributions of z.scores for each motif
par(mfrow = c(1,2))

zmat.df <- melt(as.data.frame(z.mat))
mot.pro <- ggplot(zmat.df, aes(x = variable, y = value)) + geom_boxplot()
mot.pro + geom_hline(aes(yintercept = 0))
ggsave(filename = "motif_profile1.jpeg", height = 15, width = 15, units = "cm")

zmatnorm.df <- melt(as.data.frame(z.mat.norm))
mot.pro.norm <- ggplot(zmatnorm.df, aes(x = variable, y = value), na.rm = T) + geom_boxplot()
mot.pro.norm + geom_hline(aes(yintercept = 0))
ggsave(filename = "motif_profile2.jpeg", height = 15, width = 15, units = "cm")

###########
## Allesina and Pascual 2008 data
apwebs <- c("benguela", "bridgebrook", "broom", "chesapeake", "littlerock", "reef", "shelf", "stmarks", "stmartin")
ap.motifs <- z.mat[c(3, 6, 8, 13, 31, 38, 39, 41,42),]
rownames(ap.motifs) <- apwebs
ap.qss <- matrix(c(.79, .86, .96, .97, .80, .75, .72, .89, .87), ncol = 1)
ap.motif.qss <- as.data.frame(cbind(ap.motifs, ap.qss))

fittest <- lm(ap.qss ~ s1 + s2 + s3 + s4 + s5, data = ap.motif.qss)
summary(fittest)

nulldistr.web <- do.call(rbind, websplit)


###  SAVE WORKSPACE 
setwd("~/Desktop/GitHub/Ecological Networks/Projects/Rdata")
save.image("foodweb.Rdata")