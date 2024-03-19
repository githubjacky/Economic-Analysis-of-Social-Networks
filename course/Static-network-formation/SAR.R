rm(list=ls())
source("C:/Users/ssunr/Dropbox/2022_NCKU_camp/code/SCSAR/DGP.r")
library(mvtnorm)
library(tictoc)
rm(list=ls()[ls()!="X" & ls()!="Y" & ls()!="W" & ls()!="C"])


N <- 30
G <- 50
T <- 3000    # number of iterations during Markov process
R <- 10       # number of replication for Monte Carlo experiment

lambda_R <- t(rep(0, R))
beta_R <- matrix(0, nrow = 2, ncol = R)
sige_R <- t(rep(0, R))

for (r in 1:R) {
  print(paste0('r', r), quote = FALSE)
  
  start_time <- Sys.time()
  
  beta_0 <- c(0,0)
  B_0 <- diag(2)*3
  inv_B_0 <- solve(B_0)
  ALPHA_0 <- 1
  sig0 <- c(1,0.5)
  rho_0 <- 2.2
  eta_0 <- 0.1
  
  c_1 <- 1
  c_2 <- 0.1
  c_3 <- 0.1
  c_4 <- 1
  
  acc_1 <- 0
  acc_2 <- 0
  acc_3 <- 0
  acc_4 <- rep(0,G)
  
  acc_rate1 <- matrix(0, nrow = T, ncol = 1)
  acc_rate2 <- matrix(0, nrow = T, ncol = 1)
  acc_rate3 <- matrix(0, nrow = T, ncol = 1)
  acc_rate4 <- matrix(0, nrow = G, ncol = T)
  
  lambda_T <- t(rep(0,T))                     # save for lambda
  beta_T   <- matrix(0, nrow = 2, ncol = T)   # save for betas
  alpha_T  <- matrix(0, nrow = G, ncol = T)   # save for alpha
  sige_T   <- t(rep(0,T))                     # save for Sigma_e^2
  
  # starting value of draw
  sige_T[1] = 1
  
  for (t in 2:T) {
    accept = 0
    # propose lambda^*
    if(t <= 2){
      lambda_1 = rnorm(1,lambda_T[t-1],0.1)
    } else{
      lambda_1 = rnorm(1, lambda_T[t-1], sd(lambda_T[1:(t-1)])*2.38)*0.95 + rnorm(1, lambda_T[t-1], 0.1)*0.05
    }
    
    
    pp_l <- 1
    V <- sige_T[t-1]*diag(N)
    inv_V <- solve(V)
    
    for (g in 1:G) {
      S_1 <- diag(N)-lambda_1*W[[r]][,,g]
      S_2 <- diag(N)-lambda_T[t-1]*W[[r]][,,g]
      ep_1 <- S_1 %*% Y[[r]][,g] - (cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t-1]-rep(1,N) * alpha_T[g,t-1])
      ep_2 <- S_2 %*% Y[[r]][,g] - (cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t-1]-rep(1,N) * alpha_T[g,t-1])
      like_1 <- det(S_1)*exp(-(1/2)*(t(ep_1) %*% inv_V %*% ep_1)) 
      like_2 <- det(S_2)*exp(-(1/2)*(t(ep_2) %*% inv_V %*% ep_2)) 
      pp_l <- pp_l*(like_1/like_2)      
    }
    pp_l = min(pp_l,1)
    
    
    if(runif(1) <= pp_l){
      lambda_T[t] = lambda_1
      acc_2 = acc_2 + 1
    }else{
      lambda_T[t] = lambda_T[t-1]
    }
    acc_rate2[t,1] = acc_2/t
    
    # THE SAMPLING OF BETA FROM PROSTERIOR DISTRIBUTION #
    ZVY=0
    ZVX=0
    
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      YY <- SS %*% Y[[r]][,g] - rep(1,N)*alpha_T[g,t-1]
      ZZ <- cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g])
      ZVX <- ZVX + t(ZZ) %*% inv_V %*% ZZ
      ZVY <- ZVY + t(ZZ) %*% inv_V %*% YY
    }
    beta_T[,t] <- rmvnorm(1, solve(inv_B_0 + ZVX) %*% (inv_B_0 %*% beta_0 + ZVY), solve(inv_B_0 + ZVX))
    
    # THE SAMPLING OF SIGMA_E^2 FROM PROSTERIOR DISTRIBUTION #
    ep_v <- rep(0, N*G)
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      ep <- SS %*% Y[[r]][,g] - cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t] - rep(1,N)*alpha_T[g,t-1]
      ep_v[((g-1)*N+1) : (g*N)] <- ep
    }
    rho_1 <- rho_0 + length(ep_v)
    sige_T[t] <- (t(ep_v) %*% ep_v + eta_0)/rchisq(1, df = rho_1)
    
    # THE SAMPLING OF ALPHA_G FROM PROSTERIOR DISTRIBUTION #
    dd <- (ALPHA_0^(-1) + sige_T[t]^(-1)* rep(1,N) %*% diag(N) %*% rep(1,N))^(-1)
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      YY <- SS %*% Y[[r]][,g]
      XX <- cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g])
      alpha_T[g,t] <- sige_T[t]^(-1)*dd* ( rep(1,N) %*% diag(N) %*% (YY - XX %*% beta_T[,t]) ) + rnorm(1)*sqrt(dd)
    }
    
    if(((t/100)-round(t/100)) == 0){
      print(paste0('iteration: ', t), quote = FALSE)
      print(paste0(c("draw of lambda: ", lambda_T[t]), collapse=" "))
      print(paste0(c("draw of beta: ", t(beta_T[,t])), collapse=" "))
      print(paste0(c("draw of sigma: ", sige_T[t]), collapse=" "))
      
      par(mfrow = c(3,1))
      plot(lambda_T[1:t], type='l', main = expression(lambda))
      plot(beta_T[1,1:t], type='l', main = expression(beta[1]))
      plot(beta_T[2,1:t], type='l', main = expression(beta[2]))
      #dev.control(displaylist = 'enable') #p36
    }
  }
  print(paste0("lambda", mean(lambda_T[1000:t])))
  print(paste0("beta", rowMeans(beta_T[,1000:t])))
  print(paste0("sigma_e^2", mean(sige_T[1000:t])))
  
  end_time <- Sys.time()
  time_needed <- end_time - start_time
  print(paste("Round", r, "was finished after", time_needed, "mins."))
  
  lambda_R[r] = mean(lambda_T[1000:t])
  beta_R[,r] = rowMeans(beta_T[,1000:t])
  sige_R[r] = mean(sige_T[1000:t])
  
  if(r < R){
    print(paste("Press Enter to continue round", r+1))
  }else{
    print(paste("This is the last round of Monte Carlo experiment"))
  }
  
}

print(c(paste0("lambda:", apply(lambda_R,1,mean)), apply(lambda_R,1,sd))) # or rowMeans(lambda_R)
print(c(paste0("beta:", apply(beta_R,1,mean)), apply(beta_R,1,sd)))
print(c(paste0("sigma_e^2:", apply(sige_R,1,mean)), apply(sige_R,1,sd)))
