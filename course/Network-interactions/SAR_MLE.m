%*******************************************************
% ESTIMATE SAR MODEL USING MLE ESTIMATION
%*******************************************************

addpath('C:\Users\ssunr\Dropbox\teaching_NTU\Econ7217\code')
ssd = 20200322;
rng(ssd);

clear;

load('dataset_network_interactions.mat')

YY=variable(:,1);
XX=variable(:,6:17);
WW=network;

tic;
init_lambda=0.03;
init_beta=zeros(25,1);
init_sige=0.5;

init=[init_lambda; init_beta; init_sige];
options = optimoptions('fminunc','Algorithm','quasi-newton', ...
    'MaxIter',5e10,'MaxFunEvals',5e10, 'Tolfun', 1e-7, ...
    'Display','iter');


[est_SAR,like_SAR,exitflag,output,grad,hess] = fminunc(@(init)SAR_likeli(init,YY,XX,WW),init,options);


lambda=est_SAR(1);
beta=est_SAR(2:25);
sige=est_SAR(26);

%var=inv(grad*grad');
var=inv(hess);
se=diag(sqrt(var));

disp([est_SAR se]);


disp('lambda');disp([est_SAR(1) se(1)]);
disp('beta');disp([est_SAR(2:26) se(2:26)]);
disp('sigma_e^2');disp([est_SAR(27) se(27)]);



