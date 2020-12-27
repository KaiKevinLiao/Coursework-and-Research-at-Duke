cd 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Sorting\code'
load 'WI_sort_100.mat'

% Number of Maximum Iteration
options = optimset('MaxFunEvals',50, 'Display','iter'); 
% Number of Sites
numP = 100;
% Number of Individuals
numI = 2404;
% Number of Monto-Carlos Simulation
numK = 50;

%%%%%%%%%%%%%%
% Initialization of Parameters
%%%%%%%%%%%%%%

alpha_A = -0.1; 
gamma_A = [0,0.4,1,0.4];
sigma_v = 1;

para = [alpha_A, gamma_A,sigma_v]';

%%%%%%%%%%%%%%
% Estimation
%%%%%%%%%%%%%%

f = @(para)LiFunc_6(para,numP,numI,numK,L,Z,X,I);

% Estimation through MLE
[para_star, fval, history] = fminsearch(f,para,options);


%%%%%%%%%%%%%%
% Recover theta using alpha_hat and gamma_hat
%%%%%%%%%%%%%%


% Read Results
alpha_hat = para_star(1) 
gamma_hat = para_star(2:5) 
sigma_v_hat = para_star(6)

sigma = X(6,:);            
s = zeros(1,numP);         
expm = zeros(numI,numP);    
ind_llh = zeros(numI,numP); 

exp_temp = zeros(numK,numP);
llh_temp = zeros(numK,numP);
denom_temp = zeros(numK,1);
v = randn(numI,numK) * sigma_v_hat;   

theta_A = zeros(numP,1);
theta = ones(numP,1);

 while sum((theta - theta_A).^2) > 0.00001  
        
        theta = theta_A; % update theta
        for i = 1:numI  
            for k = 1:numK
                for j = 1:numP
                    exp_temp(k,j) = exp(theta(j,1) + alpha_hat * Z(i,j) +...,
                        gamma_hat(1) * X(5,j) * I(2,i)+..., 
                        gamma_hat(2) * X(2,j) * I(2,i)+..., 
                        gamma_hat(3) * X(1,j) * I(1,i)+..., 
                        gamma_hat(4) * X(3,j) * I(1,i)+..., 
                        v(i,k)   *  X(3,j));    
                end
            end
            denom = sum(exp_temp,2);
            for k = 1:numK      
                for j = 1:numP
                    llh_temp(k,j) = exp_temp(k,j) / denom(k,1);
                end
            end
            ind_llh(i,:) = mean(llh_temp);
        end
        
        s = sum(ind_llh) / numI;        
        
        % Update theta
        
        for j = 1:numP
            theta_A(j,1) = theta(j,1) + (log(sigma(1,j)) - log(s(1,j)));           
        end
        
    end
    % Normalization of Theta
    mean_theta = mean(theta_A);
    theta_A = theta_A - mean_theta; 

X_1=[ones(100,1),X(1:5,:)'];  
beta = regress(theta_A,X_1)
