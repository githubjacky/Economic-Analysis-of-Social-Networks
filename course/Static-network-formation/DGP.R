rm(list=ls())
setwd("C:/Users/ssunr/Dropbox/2022_NCKU_camp/code/SCSAR")
#install.packages("tictoc")
#install.packages("mvtnorm")
library(mvtnorm)
library(tictoc)
tic("Time")
set.seed(20221117)

# ******* DGP Paremeters ******* #
gamma <- c(-1.5, 0.5, 1)
beta <- c(0.5, 0.5)
Sigma <- matrix(c(1, 0.75, 0.75, 1.25), ncol = 2)
sigma_a <- 0.05
lambda <- 0.05
delta <- 0.3
mu <- c(0,0)

# choose group size and number of groups #
N <- 30       
G <- 50
R <- 10

C <- array(list(), R)
W <- array(list(), R)
X <- array(list(), R)
Y <- array(list(), R)


for (r in 1:R) {

  print(r)

  z <- matrix(0, nrow = N, ncol = G)
  e <- matrix(0, nrow = N, ncol = G)
  c_1 <- matrix(runif(N*G), nrow = N, ncol = G)
  c1 <- array(runif(N*N*G), dim = c(N,N,G))
  
  for (g in 1:G) {
    Mat <- mvtnorm::rmvnorm(N,mu,Sigma) # use rmvnorm in package mvtnorm
    for (i in 1:N) {
      for (j in 1:N) {
        if (c_1[i,g] >= 0.7 & c_1[j,g] >= 0.7){
          c1[i,j,g] = 1
        } else if (c_1[i,g] <= 0.3 & c_1[j,g] <= 0.3){
          c1[i,j,g] = 1
        } else{
          c1[i,j,g] = 0
        }
        
        if (i == j){
          c1[i,j,g] = 0
        }
      }
      z[i,g] <- Mat[i,1]
      e[i,g] <- Mat[i,2]
    }
  }
  
  C[[r]] <- c1 
  p <- array(0, dim = c(N,N,G))
  w <- array(0, dim = c(N,N,G))
  for (g in 1:G) {
    for (i in 1:N) {
      for (j in 1:N) {
        psi <- gamma[1] + gamma[2]*c1[i,j,g] - gamma[3]*abs(z[i,g]-z[j,g])
        p[i,j,g] <- exp(psi)/(1+exp(psi))   # Logit setting
        if (runif(1) <= p[i,j,g]){
          w[i,j,g] = 1
        } else{
          w[i,j,g] = 0
        }
        if (i == j){
          w[i,j,g] =0
        }
      }
    }
  }
  W[[r]] <- w
  x <- matrix(rnorm(N*G)*2, nrow = N, ncol = G)
  alpha <- rnorm(G)*sqrt(sigma_a) #apply(x,2,mean)*delta +
  y <- matrix(0, nrow = N, ncol = G)
  for (g in 1:G) {
    S <- diag(N) - lambda*w[,,g]
    y[,g] <- solve(S) %*% (cbind(x[,g], w[,,g] %*% x[,g]) %*% beta + rep(1,N) * alpha[g] + e[,g])
  }
  X[[r]] <- x
  Y[[r]] <- y
}


#C_file <- tempfile("C", fileext = ".rds")
#saveRDS(C, C_file, compress = FALSE)
#W_file <- tempfile("W", fileext = ".rds")
#saveRDS(W, W_file, compress = FALSE)
#X_file <- tempfile("X", fileext = ".rds")
#saveRDS(X, X_file, compress = FALSE)
#Y_file <- tempfile("Y", fileext = ".rds")
#saveRDS(Y, Y_file, compress = FALSE)

toc()

