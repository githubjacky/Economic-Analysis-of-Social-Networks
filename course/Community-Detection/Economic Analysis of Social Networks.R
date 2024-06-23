# Example and code provided by Jordi Casas-Roma 
# accompanying with paper ``Community-preserving anonymization of graphs'' 

install.packages('igraphdata')
#install.packages('igraph')

library(igraph)
require(igraphdata)
data(karate) 


CorenessLayout <- function(g) {
coreness <- graph.coreness(g);
xy <- array(NA, dim=c(length(coreness), 2));
 
shells <- sort(unique(coreness));
for(shell in shells) {
v <- 1 - ((shell-1) / max(shells));
nodes_in_shell <- sum(coreness==shell);
angles <- seq(0,360,(360/nodes_in_shell));
angles <- angles[-length(angles)]; # remove last element
xy[coreness==shell, 1] <- sin(angles) * v;
xy[coreness==shell, 2] <- cos(angles) * v;
}
return(xy);
}
 
# compute coreness of karate network
coreness <- graph.coreness(karate);
# assign colors
colbar <- rainbow(max(coreness));
# create layout
ll <- CorenessLayout(karate);
# plot
plot(karate, layout=ll, vertex.size=15, vertex.color=colbar[coreness],vertex.frame.color=colbar[coreness], main='Coreness')