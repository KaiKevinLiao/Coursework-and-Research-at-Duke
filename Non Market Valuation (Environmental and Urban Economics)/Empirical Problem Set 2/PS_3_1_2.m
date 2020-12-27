cd 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Sorting\code';
load 'WI_sort_100.mat'

% Number of Maximum Iteration
options = optimset('MaxFunEvals',3000); 
% Number of Sites
numP = 100;
% Number of Individuals
numI = 2404;

% Guess the initial value of parameters
alpha_init = -0.1;
beta_init = [0,0,0,0,0,50];
gamma_init = [0.5,0.4,0.5,0.6];
para = [alpha_init,beta_init,gamma_init]';

% Read the likelihood function
f = @(para)LiFunc_5(para,numP,numI,L,Z,X,I);

% Estimation through MLE
para_star = fminsearch(f,para,options);

% Read Results
alpha_hat = para_star(1) 
beta_hat = para_star(2:7) 
gamma_hat = para_star(8:11) 