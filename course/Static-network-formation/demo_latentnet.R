#ECON7217 Statistic Network Formation: Latent Space Model with 'latentnet'

setwd("C:/Users/ssunr/Dropbox/teaching_NTU/Econ7217/code");
install.packages("latentnet");
install.packages("ergm");
install.packages("network");
install.packages("sna");
install.packages("blockmodels")
library(igraph);
library(ergm);
library(sna);
library(latentnet);
library(network);
library(blockmodels)

load("dataset_latentnet.RData");

plot(net)
plot(net, vertex.col = "male")

#Block Model with 'sna'

eq <- equiv.clust(net);
block <- blockmodel(net, eq, k=3);
plot(block);



## Sample stochastic block model with 'igraph'
pm <- cbind( c(.1, .001, 0.03), c(0.001, 0.2, 0.05), c(0.03, 0.05, .3) )
g <- sample_sbm(100, pref.matrix=pm, block.sizes=c(20,40,40))


## estimate stochastic block model with 'blockmodels'
adj_g <-  as.matrix(as_adjacency_matrix(g))
my_model <- BM_bernoulli("SBM",adj_g)
my_model$estimate()
my_model$model_parameters
my_model$plot_obs_pred(3)
which.max(my_model$ICL)


# Latent Space Effects with "latentnet"
#
# euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL, mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL) 
# Adds a term to the model equal to the negative Eucledean distance -dist(Z[i],Z[j]), where Z[i] and Z[j] are the positions of their respective actors in an unobserved social space. 
#
# bilinear(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL, mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
# Adds a term to the model equal to the inner product of the latent positions: sum(Z[i]*Z[j]), where Z[i] and Z[j] are the positions of their respective actors in an unobserved social space.
#
# Actor-specific effects
# rsender(var=1, var.df=3)
# Random sender effect. Adds a random sender effect to the model


model1_fomula <- formula(net ~ euclidean(d=2));
model1 <- ergmm(model1_fomula,
                control=ergmm.control(burnin=100000,sample.size= 10000,interval=5));
summary(model1);
plot(model1, pie = TRUE, vertex.cex = 2.5);
#mcmc.diagnostics(model2);

model2_fomula <- formula(net ~ nodematch("white") + nodematch("male") +
                               euclidean(d=2));
model2 <- ergmm(model2_fomula,
                control=ergmm.control(burnin=100000,sample.size= 10000,interval=5));
summary(model2);
plot(model2, pie = TRUE, vertex.cex = 2.5);
#mcmc.diagnostics(model2);


