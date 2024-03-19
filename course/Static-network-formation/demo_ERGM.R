install.packages('statnet')
library(statnet)
data(package='ergm')

# tutorial using ERGM in R https://sites.duke.edu/dnac/an-ergm-tutorial-using-r/

## Example 1

data(florentine)
flomarriage
plot(flomarriage, main="Florentine Marriage", cex.main=0.8)
summary(flomarriage ~ edges+triangle)

model1 <-ergm(flomarriage ~ edges+triangle)
summary(model1)

## Example 2

data(faux.desert.high)
faux.desert.high
plot(faux.desert.high, main="Faux desert High School", cex.main=0.8)

summary(faux.desert.high ~ edges+mutual+triangle);

model1 <-ergm(faux.desert.high ~ edges + mutual + nodematch("race") + nodematch("sex"))
summary(model1)

model2 <-ergm(faux.desert.high ~ edges+mutual+triangle, control = control.ergm(seed=1, MCMC.samplesize=1000),
              verbose=T);
summary(model2)

model3 <-ergm(faux.desert.high ~ edges+mutual+triangle, control = control.ergm(seed=1, MCMC.samplesize=20000),
              verbose=T);
summary(model3)


# https://eehh-stanford.github.io/SNA-workshop/ergm-predictions.html#the-gwesp-statistic

model4 <-ergm(faux.desert.high ~ edges+mutual+gwesp(0.25,fixed=T), control=control.ergm(MCMLE.maxit= 10,MCMC.interval=2000,MCMC.samplesize=2000),
              verbose=T);

summary(model4)
mcmc.diagnostics(model4)





