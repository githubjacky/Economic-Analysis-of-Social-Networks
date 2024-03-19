function [likeli] = SAR_likeli(parm,Y,X,W)
% PURPOSE: Generate log-likelihood function for the SAR model
% with group fixed effects
% SAR Model: Y_g=lambda*W_g*Y_g+X_g*beta_1+W_g*X_g*beta_2+ell_N*alpha_g+e_g
% -----------------------------------------------------------------------------

lambda=parm(1);
beta=parm(2:26);
sige=parm(27);


likeli=0;
  
N=length(W);
V=sige*eye(N);
S=eye(N)-lambda*W;
ep=S*Y-[ones(N,1) X W*X]*beta;
lnL_g=(-N/2)*log(2*pi)-(1/2)*log(det(V))+log(det(S))-(1/2)*ep'/V*ep;
likeli=likeli-lnL_g;









