cd 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Sorting\code'
load 'WI_sort_100.mat'

% Number of Maximum Iteration
options = optimset('MaxFunEvals',500, 'Display','iter');
% Number of Sites
numP = 100;
% Number of Individuals
numI = 2404;

%%%%%%%%%%%%%%
% Initialization of Parameters
%%%%%%%%%%%%%%

alpha_A = -0.1; 
gamma_A = [0,0.4,1,0.4];

para = [alpha_A, gamma_A]';


%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%

f = @(para)LiFunc_3(para,numP,numI,L,Z,X,I);

% Estimation through MLE
[para_star, fval, history] = fminsearch(f,para,options);

% Read Results
alpha_hat = para_star(1) 
gamma_hat = para_star(2:5) 


%%%%%%%%%%%%%%
% Recover theta using alpha_hat and gamma_hat
%%%%%%%%%%%%%%

theta_A = zeros(100,1);
theta = ones(100,1);

sigma = X(6,:);            
s = zeros(1,100);          
expm = zeros(2404,100);    
denom = zeros(2404,1);      
ind_llh = zeros(2404,100); 

while sum((theta - theta_A).^2) > 0.00001  %while theta not converge

    theta = theta_A; % update theta
    for i = 1:numI         
        for j = 1:numP
            expm(i,j) = exp(theta(j,1) + alpha_hat * Z(i,j) +...,
                gamma_hat(1) * X(5,j) * I(2,i)+..., 
                gamma_hat(2) * X(2,j) * I(2,i)+..., 
                gamma_hat(3) * X(1,j) * I(1,i)+..., 
                gamma_hat(4) * X(3,j) * I(1,i));    
        end
    end
    for i = 1:numI
        denom(i,1) = sum(expm(i,:));
        for j = 1:numP
            ind_llh(i,j) = expm(i,j) / denom(i,1); 
        end
    end
    s = sum(ind_llh) / numI;        
    for j = 1:numP
        theta_A(j,1) = theta(j,1) + (log(sigma(1,j)) - log(s(1,j)));        
    end
end
mean_theta = mean(theta_A);
theta_A = theta_A - mean_theta; 


%%%%%%%%%%%
% IV Strategy
%%%%%%%%%%%

X_1=[ones(100,1),X'];  %exclude shares
beta_A = regress(theta_A,X_1);
theta = theta_A;
beta = beta_A(2:6);

exp_iv = zeros(2404,100);  
denom_iv = zeros(2404,1);     
share_iv_temp = zeros(2404,100);  

for i = 1:numI         
    for j = 1:numP
        exp_iv(i,j) = exp(beta' * X(1:5,j) + alpha_hat * Z(i,j) +...,
            gamma_hat(1) * X(5,j) * I(2,i)+..., 
            gamma_hat(2) * X(2,j) * I(2,i)+..., 
            gamma_hat(3) * X(1,j) * I(1,i)+..., 
            gamma_hat(4) * X(3,j) * I(1,i));    
    end
end
for i = 1:numI
    denom_iv(i,1) = sum(exp_iv(i,:));
    for j = 1:numP
        share_iv_temp(i,j) = exp_iv(i,j) / denom_iv(i,1); 
    end
end
share_iv = sum(share_iv_temp) / numI       

X_iv = [ones(100,1),X(1:5,:)',share_iv'];
temp = fitlm(X_iv,sigma');
sigma_iv_percent = predict(temp,X_iv)*100;


X_2 = [ones(100,1),X(1:5,:)',sigma_iv_percent];
beta_iv = regress(theta,X_2)
