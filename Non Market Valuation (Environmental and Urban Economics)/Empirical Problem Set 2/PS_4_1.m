%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation of Partial Equilibrium
%
% Though this is a partial equilibrium analysis, I still allow so kind of
% re-optimization of fishers. The fishers can resort base on the exougenous
% change. However, I do not allow "feedback effect" - the effect of the new
% congestion levels. 
%%%%%%%%%%%%%%%%%%%%%%%%%%



cd 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Sorting\code'
load 'Parameters_For_Welfare_Analysis.mat'
load 'WI_sort_100.mat'

% Number of Maximum Iteration
options = optimset('MaxFunEvals',500, 'Display','iter'); 
% Number of Sites
numP = 100;
% Number of Individuals
numI = 2404;

exp_base = zeros(numI,numP);
exp_1 = zeros(numI,numP);
exp_2 = zeros(numI,numP);
exp_3 = zeros(numI,82);
exp_4 = zeros(numI,numP);

Phi_base = [];

Phi_1 = [];

Phi_2 = [];
Phi_base_a = [];
Phi_base_ua = [];
Phi_2_a = [];
Phi_2_ua = [];

Phi_base_3 = [];
Phi_4 = [];

% Rearrage X matrix so that it is better for regression
X_temp = X;
X = [ones(numP,1),X_temp'];

% Take the unobserved term for theta
beta_iv(7) = 100 * beta_iv(7);
xi = theta - X*beta_iv;

alpha = alpha_hat;
gamma = gamma_hat;

%%%%%%%%%
% Baseline Results
%%%%%%%%%

for i = 1:numI         
    for j = 1:numP
        exp_base(i,j) = exp(theta(j,1) + alpha * Z(i,j) +...,
            gamma(1) * X_temp(5,j) * I(2,i)+..., 
            gamma(2) * X_temp(2,j) * I(2,i)+..., 
            gamma(3) * X_temp(1,j) * I(1,i)+..., 
            gamma(4) * X_temp(3,j) * I(1,i));    
    end
end
Phi_base = log(sum(exp_base')');

%%%%%%%%%
% Scenario 1
%%%%%%%%%

Walleye_temp = 1.3 * X(:,4);
X_1 = [X(:,1:3),Walleye_temp,X(:,5:7)];
theta_1 = X_1 * beta_iv + xi ;
for i = 1:numI         
    for j = 1:numP
        exp_1(i,j) = exp(theta_1(j,1) + alpha * Z(i,j) +...,
            gamma(1) * X_1(j,6) * I(2,i)+..., 
            gamma(2) * X_1(j,3) * I(2,i)+..., 
            gamma(3) * X_1(j,2) * I(1,i)+..., 
            gamma(4) * X_1(j,4) * I(1,i));    
    end
end
Phi_1 = log(sum(exp_1')');
CV_1 = (-1)/alpha * (Phi_1 - Phi_base);
CV_mean_1 = mean(CV_1)

%%%%%%%%%
% Scenario 2
%%%%%%%%%

% Construct Counterfactual X_2 and theta
X_2 = X;
for i = 1:numP
    if X(i,7) > 0.015
        X_2(i,4) = 1.3 * X(i,4);
    end
end
theta_2 = X_2 * beta_iv + xi;


% Construct the Share Matrix
for i = 1:numI         
    for j = 1:numP
        exp_2(i,j) = exp(theta_2(j,1) + alpha * Z(i,j) +...,
            gamma(1) * X_2(j,6) * I(2,i)+..., 
            gamma(2) * X_2(j,3) * I(2,i)+..., 
            gamma(3) * X_2(j,2) * I(1,i)+..., 
            gamma(4) * X_2(j,4) * I(1,i));    
    end
end

% Overall welfare effect
Phi_2 = log(sum(exp_2')');
CV_2 = (1)/alpha * (Phi_base - Phi_2);
CV_mean_2 = mean(CV_2)

% Welfare effects for affacted and unaffacted
for i = 1:numI
    C = L(i,1);
    if X_1(C,7) > 0.015
        Phi_2_a = [Phi_2_a;Phi_2(i,1)];
        Phi_base_a = [Phi_base_a;Phi_base(i,1)];
    else
        Phi_2_ua = [Phi_2_ua;Phi_2(i,1)];
        Phi_base_ua = [Phi_base_ua;Phi_base(i,1)];
    end
end

CV_2_a = (-1)/alpha * (Phi_2_a - Phi_base_a) ;
CV_mean_2_a = mean(CV_2_a)
CV_2_ua = (-1)/alpha * (Phi_2_ua - Phi_base_ua) ;
CV_mean_2_ua = mean(CV_2_ua)

%%%%%%%%%
% Scenario 3
%%%%%%%%%

X_3 = [];
xi_3 = [];
Z_3 =[];
Phi_base_3 = Phi_base_ua;

% Construct Counterfactual X_2 and theta
for i = 1:numP
    if X_1(i,7) <= 0.015
        X_3 = [X_3 ; X(i,:)];
        xi_3 = [xi_3 ; xi(i,1)];
        Z_3 = [Z_3 , Z(:,i)];
    end
end
theta_3 = X_3 * beta_iv + xi_3;

% Construct the Share Matrix
for i = 1:numI
    C = L(i,1);
    if X_1(C,7) <= 0.015
        for j = 1:82
            exp_3(i,j) = exp(theta_3(j,1) + alpha * Z_3(i,j) +...,
                gamma(1) * X_3(j,6) * I(2,i)+..., 
                gamma(2) * X_3(j,3) * I(2,i)+..., 
                gamma(3) * X_3(j,2) * I(1,i)+..., 
                gamma(4) * X_3(j,4) * I(1,i));    
        end
    end
end

% Welfare effects
rowsToDelete = all(exp_3==0, 2);
exp_3(rowsToDelete,:)=[];
Phi_3 = log(sum(exp_3')');
CV_3 = (-1)/alpha * (Phi_3 - Phi_base_ua) ;
Mean_CV_3 = mean(CV_3)


%%%%%%%%%
% Scenario 4
%%%%%%%%%

Z_4 = Z;
for i = 1:numP
    if X(i,7) > 0.015
        Z_4(:,i) = 10 + Z_4(:,i);
    end
end
for i = 1:numI         
    for j = 1:numP
        exp_4(i,j) = exp(theta(j,1) + alpha * Z_4(i,j) +...,
              gamma(1) * X_temp(5,j) * I(2,i)+..., 
              gamma(2) * X_temp(2,j) * I(2,i)+..., 
              gamma(3) * X_temp(1,j) * I(1,i)+..., 
              gamma(4) * X_temp(3,j) * I(1,i));    
    end
end
Phi_4 = log(sum(exp_4')');
Phi_4_a = [];
Phi_4_ua = [];
for i = 1:numI
    C = L(i,1);
    if X_1(C,7) > 0.015
        Phi_4_a = [Phi_4_a;Phi_4(i,1)];
    else
        Phi_4_ua = [Phi_4_ua;Phi_4(i,1)];
    end
end
CV_4_a = (-1)/alpha * (Phi_4_a - Phi_base_a) ;
CV_4_ua = (-1)/alpha * (Phi_4_ua - Phi_base_ua) ;
Mean_CV_4_a = mean(CV_4_a)
Mean_CV_4_ua = mean(CV_4_ua)





