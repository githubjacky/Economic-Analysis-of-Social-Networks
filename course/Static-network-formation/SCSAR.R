#*******************************************************
# ESTIMATE SC-SAR MODEL USING BAYESIAN ESTIMATION
#*******************************************************

source("C:/Users/ssunr/Dropbox/2022_NCKU_camp/code/SCSAR/DGP.r")
library(mvtnorm)
rm(list=ls()[ls()!="X" & ls()!="Y" & ls()!="W" & ls()!="C"])

set.seed(20200322)

N=30
G=50
T=11000   # number of discarded samples during Markov process
R=10      # number of replication for Monte Carlo experiment

lambda_R <- t(rep(0, R))
beta_R <- matrix(0, nrow = 2, ncol = R)
sige_R <- t(rep(0, R))
sigez_R <- t(rep(0, R))
gamma_R <- matrix(0, nrow = 3, ncol = R)

for (r in 1:R) {
  print(paste0('r', r), quote = FALSE)
  
  start_time <- Sys.time()
  
  c_1 <- 1
  c_2 <- 0.1
  c_3 <- 0.1
  c_4 <- 1
  
  acc_1 <- 0
  acc_2 <- 0
  acc_3 <- 0
  acc_4 <- rep(0,G)
  
  # assign parameter in prior distributions #
  gamma_0 <- c(0,0,0)
  beta_0 <- c(0,0)
  G_0 <- diag(3)*3
  B_0 <- diag(2)*3
  inv_B_0 <- solve(B_0)
  rho_0 <- 2.2
  eta_0 <- 0.1
  sig0 <- c(1,0.5)
  ALPHA_0 <- 100
  
  acc_rate1 <- matrix(0, nrow = T, ncol = 1)
  acc_rate2 <- matrix(0, nrow = T, ncol = 1)
  acc_rate3 <- matrix(0, nrow = T, ncol = 1)
  acc_rate4 <- matrix(0, nrow = G, ncol = T)

  zz_T <- matrix(0, nrow = N, ncol = G)       # save for Z
  gamma_T <- matrix(0, nrow = 3, ncol = T)    # save for gamma
  lambda_T <- t(rep(0,T))                     # save for lambda
  beta_T   <- matrix(0, nrow = 2, ncol = T)   # save for betas
  alpha_T  <- matrix(0, nrow = G, ncol = T)   # save for alpha
  sige_T   <- t(rep(0,T))                     # save for Sigma_e^2
  sigez_T   <- t(rep(0,T))                    # save for Sigma_ez
  
  # starting value of draw
  gamma_T[,1] = c(-1.3,0.5,1)
  sige_T[1] = 1
  sigez_T[1] = 0.5
  
  for (t in 2:T) {
    accept = 0
    # propose lambda^*
    if(t < 100){
      lambda_1 = rnorm(1,lambda_T[t-1],0.1)
    } else{
      lambda_1 = rnorm(1, lambda_T[t-1], sd(lambda_T[1:(t-1)])*2.38)*0.95 + rnorm(1, lambda_T[t-1], 0.1)*0.05
    }
   
    
    if (t < 100){
      gamma_1 = rmvnorm(1, t(gamma_T[,t-1]), diag(3)*0.1^2/3)
    }else{
      gamma_1 = rmvnorm(1, t(gamma_T[,t-1]), cov(t(gamma_T[,1:(t-1)]))*2.38^2/3)*0.95 + rmvnorm(1,t(gamma_T[,t-1]), diag(3)*0.1^2/3)*0.05
    }
    
    pp_G = 1
    pp_l = 1
    
    # THE M-H ALGORITHM FOR SAMPLING Z #
    
    V <- (sige_T[t-1]-sigez_T[t-1]^2)*diag(N)
    inv_V <- solve(V)
    
    for (g in 1:G) {
      zz_1 <- zz_T[,g]
      acc_4v <- 0
      for (v in 1:N) {
        zz_1[v] <- rnorm(1)*c_4 + zz_T[v,g]
        pp <- 1
        for (i in 1:N) {
          if(i == v){
            for (j in 1:N) {
              if(j == v){
                p_1 = 1
                p_2 = 1
              }else{
                psi_1 <- gamma_T[1,t-1] + gamma_T[2,t-1]*C[[r]][i,j,g] - gamma_T[3,t-1]*abs(zz_1[i]-zz_1[j])
                psi_2 <- gamma_T[1,t-1] + gamma_T[2,t-1]*C[[r]][i,j,g] - gamma_T[3,t-1]*abs(zz_T[i,g]-zz_T[j,g])
                p_1 <- exp(psi_1*W[[r]][i,j,g])/(1+exp(psi_1))
                p_2 <- exp(psi_2*W[[r]][i,j,g])/(1+exp(psi_2))
              }
              pp <- pp*(p_1/p_2)
            }
          }else{
            psi_1 <- gamma_T[1,t-1] + gamma_T[2,t-1]*C[[r]][i,v,g] - gamma_T[3,t-1]*abs(zz_1[i]-zz_1[v])
            psi_2 <- gamma_T[1,t-1] + gamma_T[2,t-1]*C[[r]][i,v,g] - gamma_T[3,t-1]*abs(zz_T[i,g]-zz_T[v,g])
            p_1 <- exp(psi_1*W[[r]][i,v,g])/(1+exp(psi_1))
            p_2 <- exp(psi_2*W[[r]][i,v,g])/(1+exp(psi_2))
            pp <- pp*(p_1/p_2)
          }
        }
        SS <- diag(N) - lambda_T[t-1]*W[[r]][,,g]
        ep <- SS %*% Y[[r]][,g] - cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t-1] - rep(1,N)*alpha_T[g,t-1] 
        like_Y1 <- exp(-(1/2)*t(ep-sigez_T[t-1]*zz_1) %*% inv_V %*% (ep-sigez_T[t-1]*zz_1))
        like_Y2 <- exp(-(1/2)*t(ep-sigez_T[t-1]*zz_T[,g]) %*% inv_V %*% (ep-sigez_T[t-1]*zz_T[,g]))
        pp <- pp*(like_Y1/like_Y2)*(dmvnorm(zz_1[v])/dmvnorm(zz_T[v,g]))
        pp <- min(pp,1)
        
	  if(runif(1) <= pp){
          zz_T[v,g] = zz_1[v]
          acc_4v = acc_4v + 1
        }
        zz_1 = zz_T[,g]
      }
      if(acc_4v >= N/2){
        acc_4[g] = acc_4[g] + 1
      }
      acc_rate4[g,t] = acc_4[g]/t

      # THE M-H ALGORITHM FOR SAMPLING GAMMA AND LAMBDA #
      pp = 1
      for (i in 1:N) {
        for (j in 1:N) {
          if(i == j){
            p_1 = 1
            p_2 = 1
          }else{
            psi_1 <- gamma_1[1] + gamma_1[2]*C[[r]][i,j,g] - gamma_1[3]*abs(zz_T[i,g]-zz_T[j,g])
            psi_2 <- gamma_T[1,t-1] + gamma_T[2,t-1]*C[[r]][i,j,g] - gamma_T[3,t-1]*abs(zz_T[i,g]-zz_T[j,g])
            p_1 <- exp(psi_1*W[[r]][i,j,g])/(1+exp(psi_1))
            p_2 <- exp(psi_2*W[[r]][i,j,g])/(1+exp(psi_2))
          }
          pp <- pp*(p_1/p_2)
        }
      }
      pp_G = pp_G*pp
      S_1 <- diag(N)-lambda_1*W[[r]][,,g]
      S_2 <- diag(N)-lambda_T[t-1]*W[[r]][,,g]
      ep_1 <- S_1 %*% Y[[r]][,g] - (cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t-1]-rep(1,N) * alpha_T[g,t-1])
      ep_2 <- S_2 %*% Y[[r]][,g] - (cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t-1]-rep(1,N) * alpha_T[g,t-1])
      like_1 <- (det(S_1))*exp(-(1/2)*(t(ep_1 - sigez_T[t-1]*zz_T[,g]) %*% inv_V %*% (ep_1 - sigez_T[t-1]*zz_T[,g])))
      like_2 <- (det(S_2))*exp(-(1/2)*(t(ep_2 - sigez_T[t-1]*zz_T[,g]) %*% inv_V %*% (ep_2 - sigez_T[t-1]*zz_T[,g])))
      pp_l <- pp_l*(like_1/like_2)
    }
    pp_G = pp_G*(dmvnorm(gamma_1, gamma_0, G_0)/dmvnorm(t(gamma_T[,t-1]), gamma_0, G_0))
    pp_G = min(pp_G, 1)
    if(runif(1) <= pp_G){
      gamma_T[,t] = gamma_1
      acc_1 = acc_1 + 1
    }else{
      gamma_T[,t] = gamma_T[,t-1]
    }
    acc_rate1[t,1] = acc_1/t
    
    pp_l = min(pp_l,1)

    if(runif(1) <= pp_l){
      lambda_T[t] = lambda_1
      acc_2 = acc_2 + 1
    }else{
      lambda_T[t] = lambda_T[t-1]
    }
    acc_rate2[t,1] = acc_2/t
    if(mean(acc_rate4[,t]) < 0.4){
      c_4 = c_4/1.01
    }
    if(mean(acc_rate4[,t]) > 0.6){
      c_4 = c_4*1.01
    }
    
    # THE SAMPLING OF BETA FROM PROSTERIOR DISTRIBUTION #
    ZVY=0
    ZVX=0
    
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      YY <- SS %*% Y[[r]][,g] - sigez_T[t-1]*zz_T[,g] - rep(1,N)*alpha_T[g,t-1]
      ZZ <- cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g])
      ZVX <- ZVX + t(ZZ) %*% inv_V %*% ZZ
      ZVY <- ZVY + t(ZZ) %*% inv_V %*% YY
    }
    beta_T[,t] <- rmvnorm(1, solve(inv_B_0 + ZVX) %*% (inv_B_0 %*% beta_0 + ZVY), solve(inv_B_0 + ZVX))
    
    # THE SAMPLING OF SIGMA_E^2 & SIGMA_EZ FROM PROSTERIOR DISTRIBUTION #
    accept = 0
    while(accept == 0){
      if(t <= 4){
        sig = rmvnorm(1, c(sige_T[t-1], sigez_T[t-1]), diag(2)*0.1^2)
      }else{
        sig = rmvnorm(1, c(sige_T[t-1], sigez_T[t-1]), cov(cbind(sige_T[1:t-1], sigez_T[1:t-1]))*2.38^2/2)*0.95 + 
          rmvnorm(1, c(sige_T[t-1],sigez_T[t-1]), diag(2)*0.1^2)*0.05
      }
      sige_1 = sig[1]
      sigez_1 = sig[2]
      if(sigez_1^2 < sige_1 & sigez_1 > 0){
        accept = 1
      }
    }
    
    V1 <- (sige_1 - sigez_1^2)*diag(N)
    V2 <- (sige_T[t-1] - sigez_T[t-1]^2)*diag(N)
    pp_sig <- 1
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      ep <- SS %*% Y[[r]][,g] - cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g]) %*% beta_T[,t] - rep(1,N)*alpha_T[g,t-1] 
      like_1 <- det(V1)^-(1/2)*exp((-1/2)*t(ep - sigez_1*zz_T[,g]) %*% solve(V1) %*% (ep - sigez_1*zz_T[,g]))
      like_2 <- det(V2)^(-1/2)*exp((-1/2)*t(ep - sigez_T[t-1]*zz_T[,g]) %*% solve(V2) %*% (ep-sigez_T[t-1]*zz_T[,g]))
      pp_sig = pp_sig*(like_1/like_2)
    }
    pp_sig <- pp_sig*(dmvnorm(sig,sig0)/dmvnorm(c(sige_T[t-1],sigez_T[t-1]),sig0))
    pp_sig <- min(pp_sig,1)
    if(runif(1) <= pp_sig){
      sige_T[t] = sige_1
      sigez_T[t] = sigez_1
      acc_3 = acc_3 + 1
    }else{
      sige_T[t] = sige_T[t-1]
      sigez_T[t] = sigez_T[t-1]
    }
    acc_rate3[t,1] = acc_3/t
    
    # THE SAMPLING OF ALPHA_G FROM PROSTERIOR DISTRIBUTION #
    dd <- (ALPHA_0^(-1) + (sige_T[t] - sigez_T[t]^2)^(-1)*rep(1,N) %*% diag(N) %*% rep(1,N))^(-1)
    for (g in 1:G) {
      SS <- diag(N) - lambda_T[t]*W[[r]][,,g]
      YY <- SS %*% Y[[r]][,g] - sigez_T[t]*zz_T[,g]
      XX <- cbind(X[[r]][,g], W[[r]][,,g] %*% X[[r]][,g])
      alpha_T[g,t] <- (sige_T[t] - sigez_T[t]^2)^(-1)*dd* ( rep(1,N) %*% diag(N) %*% (YY - XX %*% beta_T[,t]) ) + rnorm(1)*sqrt(dd)
    }

    if(t/10 - round(t/10) == 0){      
      print(paste0('iteration: ', t), quote = FALSE)
      print(paste0(c("draw of gamma: ", t(gamma_T[,t])), collapse=" "))
      print(paste0(c("draw of lambda: ", lambda_T[t]), collapse=" "))
      print(paste0(c("draw of beta: ", t(beta_T[,t])), collapse=" "))
      print(paste0(c("draw of sigma: ", cbind(sige_T[t], sigez_T[t])), collapse=" "))
      print(paste0(c("acc rate4: ", cbind(acc_rate4[t], c_4)), collapse=" "))

     
      par(mfrow = c(3,2))
      plot(lambda_T[1:t], type='l', main = expression(lambda))
      plot(beta_T[1,1:t], type='l', main = expression(beta[1]))
      plot(beta_T[2,1:t], type='l', main = expression(beta[2]))
      plot(gamma_T[1,1:t], type='l', main = expression(gamma[1]))
      plot(gamma_T[2,1:t], type='l', main = expression(gamma[2]))
      plot(gamma_T[3,1:t], type='l', main = expression(gamma[3]))
      #dev.control(displaylist = 'enable') #p36
    }
  }
  
  print(paste0("gamma", rowMeans(gamma_T[,1000:t])))
  print(paste0("lambda", mean(lambda_T[1000:t])))
  print(paste0("beta", rowMeans(beta_T[,1000:t])))
  print(paste0("sigma_e^2", mean(sige_T[1000:t])))
  print(paste0("sigma_ez", mean(sigez_T[1000:t])))
  
  end_time <- Sys.time()
  time_needed <- end_time - start_time
  print(paste("Round", r, "was finished after", time_needed, "mins."))
  
  gamma_R[,r] = rowMeans(gamma_T[,1000:t])
  lambda_R[r] = mean(lambda_T[1000:t])
  beta_R[,r] = rowMeans(beta_T[,1000:t])
  sige_R[r] = mean(sige_T[1000:t])
  sigez_R[r] = mean(sigez_T[1000:t])
  
  if(r < R){
    print(paste("Press Enter to continue round", r+1))
    scan(quiet=TRUE)
  }else{
    print(paste("This is the last round of Monte Carlo experiment"))
  }
  
}

print(paste0("lambda:", c(apply(lambda_R,1,mean), apply(lambda_R,1,sd)))) # or rowMeans(lambda_R)
print(paste0("beta:", c(apply(beta_R,1,mean), apply(beta_R,1,sd))))
print(paste0("gamma:", c(apply(gamma_R,1,mean), apply(gamma_R,1,sd))))
print(paste0("sigma_e^2:", c(apply(sige_R,1,mean), apply(sige_R,1,sd))))
print(paste0("sigma_ez:", c(apply(sigez_R,1,mean), apply(sigez_R,1,sd))))