install.packages('visNetwork')
install.packages('magrittr')
install.packages('data.table')

library(igraph)
library(magrittr)
library(visNetwork)
library(data.table)


setwd("C:/Users/ssunr/Dropbox/teaching_NTU/Econ7217/code")

# Example 1  Bank transfer network 

mydata_bank_transfer <- read.csv("Community Detection and Visualization in R.csv", head = T)

# Move originator and beneficiary to the first and second columns.

move = names(mydata_bank_transfer)[names(mydata_bank_transfer) == "originator" | names(mydata_bank_transfer)== "beneficiary"]
setcolorder(mydata_bank_transfer, c(move, setdiff(names(mydata_bank_transfer), move))) 


# Simple graphs are graphs which do not contain loop and multiple edges.

graph <- graph.data.frame(mydata_bank_transfer, directed=F)
graph <- simplify(graph)  
plot(graph,vertex.label.cex=0.7)


# Girvan_Newman algorithm (divisive algorithm based on betweenness centrality)

comm_Girvan_Newman <- cluster_edge_betweenness(graph,edge.betweenness = TRUE)

V(graph)$comm_Girvan_Newman <- comm_Girvan_Newman$membership

V(graph)$color <- V(graph)$comm_Girvan_Newman + 1

plot(graph,vertex.label.cex=0.7)

# plot dendogram
plot_dendrogram(comm_Girvan_Newman)


# Louvain method

comm_Louvain <- cluster_louvain(graph)

str(comm_Louvain)

V(graph)$comm_Louvain <- comm_Louvain$membership

V(graph)$color <- V(graph)$comm_Louvain + 1

plot(graph,vertex.label.cex=0.7)



# Walktrap algorithm

comm_walktrap <- cluster_walktrap(graph)

V(graph)$comm_walktrap<- comm_walktrap$membership

V(graph)$color <- V(graph)$comm_walktrap + 1

dev.new()
plot(graph,vertex.label.cex=0.7)


# Label propagation algorithm

comm_LP <- cluster_label_prop(graph)

V(graph)$comm_LP<- comm_LP$membership

V(graph)$color <- V(graph)$comm_LP + 1

dev.new()
plot(graph,vertex.label.cex=0.7)


# Infomap algorithm

comm_infomap <- cluster_infomap(graph)
imc <- cluster_infomap(graph)
membership(imc)
communities(imc)

dev.new()
plot(comm_infomap, graph, vertex.label.cex=0.7)


# Example 2 GitHub developers social network (edges are mutual follower relationships)

setwd("C:/Users/ssunr/Dropbox/teaching_NTU/Econ7217/code/git_web_ml")

mydata_Github <- read.csv("musae_git_edges.csv", head = T)

graph <- graph.data.frame(mydata_Github, directed=F)

comm_Louvain <- cluster_louvain(graph)

V(graph)$comm_Louvain <- comm_Louvain$membership

sizes(comm_Louvain)



