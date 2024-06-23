install.packages('igraph')
library(igraph)

# for reproducibility of graphs plots (plot.igraph uses random numbers)
set.seed(123)

# create an example graph
D <- read.table(header=T,text=
'from to
A B
A C
C D
C F
C E
D E
D F
E F')

g1 <- graph.data.frame(D,directed=F)

# plot the original graph
dev.new()
plot(g1)
 

# find all the largest cliques (returns a list of vector of vertiex ids)
a <- largest.cliques(g1)

# let's just take the first of the largest cliques
# (in this case there's just one clique)
clique1 <- a[[1]]

# subset the original graph by passing the clique vertices
g2 <- induced.subgraph(graph=g1,vids=clique1)

# plot the clique
dev.new()
plot(g2)
 