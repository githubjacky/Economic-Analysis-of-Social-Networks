install.packages("igraph")
library(igraph)
setwd("/Users/jackyyeh/github/Economic-Analysis-of-Social-Networks/course/Introduction-and-characterization-of-social-networks")



# reference
# https://rstudio-pubs-static.s3.amazonaws.com/74248_3bd99f966ed94a91b36d39d8f21e3dc3.html
# http://www.mjdenny.com/Preparing_Network_Data_In_R.html
# https://rdatamining.wordpress.com/2012/05/17/an-example-of-social-network-analysis-with-r-using-package-igraph/


# Example 1. generate a network object directly from a matrix

W=matrix( c(0,1,0,0,0,1,0,1,1,0,0,1,0,1,0,0,0,1,0,0,1,1,0,0,0), nrow=5, ncol=5)
g <- graph.adjacency(W)
plot(g,edge.arrow.size = 0.1)


# Example 2. generate a network object from an edge list

edgelist <- data.frame(
  source=c("A","B", "A", "A", "A","F", "B"),
  target=c("B","A", "C", "D", "F","A","E")
)

g <- graph_from_data_frame(d=edgelist, directed=F) 

plot(g)


# Example 3. generate a network object from an edge list

g <- graph.formula(1-2, 1-3, 2-3, 2-4, 3-5, 4-5, 4-6, 4-7, 5-6, 6-7);

plot(g)

V(g)

E(g)

get.edgelist(g)


# Example 4. 

dat=read.csv("sample_adjmatrix.csv",header=TRUE,row.names=1,check.names=FALSE) # choose an adjacency matrix from a .csv file
m=as.matrix(dat) # coerces the data set as a matrix
g=graph.adjacency(m,mode="undirected",weighted=NULL) # this will create an 'igraph object'



plot(g)


W2=m%*%m        # W2 shows the number of walks between every 2 nodes with 2 steps # 
W3=m%*%m%*%m    # W3 shows the number of walks between every 2 nodes with 3 steps # 

deg <- degree(g)
bet <- betweenness(g, directed = FALSE)
clo <- closeness(g, mode = "out")
evc <- evcent(g, scale = TRUE, weights = NULL)
pag <- page_rank(g)


# Calculate the average path length of the network 
average.path.length(g, directed=FALSE, unconnected=TRUE)

# Calculate the assortativity coefficient
assortativity_degree(g, directed = FALSE)

# Calculate the diameter of the network 
diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)


#save png/jpeg/pdf of the graph

png(filename = "sample_adjmatrix.png", width = 650, height = 480);
par(mfrow = c(1,2));
igraph.options(vertex.size = 7, vertex.label = NA, edge.arrow.size = 0.02,
               vertex.color = "blue");
plot(g, layout = layout.circle);
title("Sample Network");
plot(g, layout = layout.fruchterman.reingold);
title("Sample Network");
dev.off()


install.packages("igraphdata");

library(igraphdata)
data(yeast)
par(mfrow = c(1,2))
igraph.options(vertex.size = 6, vertex.label = NA, edge.arrow.size = 0.02)
plot(yeast, layout = layout.kamada.kawai)
title("yeast network")
plot(yeast, layout = layout.drl)
title("yeast network")

degree.yeast <- degree(yeast)
hist(degree.yeast,100, col = "blue",
     xlab = "Degree", ylab = "Frequency",
     main = "Degree Distribution of Yeast Dataset")

dd.yeast <- degree.distribution(yeast)
d <- 1:max(degree.yeast)-1
ind <- (dd.yeast != 0)
plot(d[ind], dd.yeast[ind], log = "xy", col = "blue",pch = 16,
     xlab = "Degree", ylab = "Frequency",
     main = "Log Degree Distribution")
 