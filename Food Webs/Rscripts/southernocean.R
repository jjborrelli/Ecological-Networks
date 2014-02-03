s.ocean <- read.csv("http://esapubs.org/archive/ecol/E092/097/diet.csv")

colnames(s.ocean)

levels(s.ocean$LOCATION)
m <- split(s.ocean, f = s.ocean$LOCATION)



