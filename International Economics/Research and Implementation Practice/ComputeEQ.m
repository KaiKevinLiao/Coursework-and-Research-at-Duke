function ComputeEQ()
% This function computes the equilibrium at parameter
%      values given by the structure param
global gamma nu mu pi NX_Data big_Y_Data Unemp_Data Emp_Data ...
    ShareEmp_Data big_L_bar_Data w_bar_Data CoeffVar_w_Data s_tilde_Data...
    Ergodic_Dist_Data big_L_bar World_big_Y_Data NX_Data;
global COUNTRY_EST K N numP NP NP_fixed TempSimAnn FUNCTION_ITERATION;
global s s_tilde Unemp ShareEmp w_bar Var_w  labor_share big_L  lambda_tilde
global x_bar w_tilde q_theta theta G_x_bar dist_out;
global b eta kappa_tilde sigma_x chi C_US C;

global lambda delta xsi zeta beta;

global big_Y;

global toler_U maxiter_U toler_out maxiter_out step_L_ini step_x_ini dist_out_vec;

global param_in;


% maximum # of iterations -- big_U loop
toler_U   = 0.00000001;
% convergence tolerance -- big_U loop
maxiter_U   = 50000;
% maximum # of iterations -- big_L and x_bar loop
toler_out = 0.0001;
% convergence tolerance -- big_L and x_bar loop
maxiter_out = 50000;
% initial steps for updating equilibrium variables
step_L_ini = 0.01;
step_x_ini = 0.01;
% storage of the dist_out for each iteration -- big_L and x_bar loop
dist_out_vec = zeros(1000000,1);

global flag_abort;
flag_abort = 0;

% need to delete afterwards
FUNCTION_ITERATION = 1;

%%%%%%%
% read parameters from param_in
%%%%%%%

% Read b
for n = 1:N
    b(n) = param_in(n);
end

% Read kappa_tilde
for k = 1:K
    for n = 1:N
        kappa_tilde(k, n) = param_in(6+(n-1)*K+k);
    end
end

% Read sigma_x
sigma_x = param_in(43);

% Read chi
% chi(sector, country) = chi(country) + chi(sector)
chi_country = zeros(N);
for n = 1:N
    chi_country(n) = param_in(43+n);
end
% Read chi_sector
chi_sector = zeros(K);
for k = 1:K
    chi_sector(k) = param_in(49+n);
end
% construct chi
for k = 1:K
    for n = 1:N
        chi(k, n) = chi_sector(k)+chi_country(n);
    end
end


% Read eta
for k = 1:K
    for n = 1:N
        eta(k, n) = param_in(55+(n-1)*K+k);
    end
end

% Read C_US
for k = 1:K
    for l = 1:K
        C_US(k,l) = param_in(91+(k-1)*K+l);
    end
end

% Construct C for all countries
for k = 1:K
    for l = 1:K
        for n = 1:N
            C(k, l, n) = C_US(k, l);
        end
    end
end


%%%%%%%
%STEP 3
%%%%%%%

%STEP 3.1 Get the value of omegabar  
for k = 1:K
    for i = 1:N
        varpi(k, i) = ((1-(1-chi(k, i))*delta)*kappa_tilde(k, i))/(delta*(1-beta));
    end
end
varpi = ((1-(1-chi)*delta).*kappa_tilde)./(delta*(1-beta));
%STEP 3.2 Check if omega_bar satisfies q(theta) <= 0
%   Detailed explanation goes here

varpi


%%%%%%%
%STEP 4 
%%%%%%%

%   return the upperbound for x, if not, return the position where this no
%   feasiable numerical solution
no_solution = zeros(K,N);
underline_x_ub = zeros(K,N);
parfor k = 1:K
    for i = 1:N
        x = sym('x');
        varpi_temp = varpi(k, i);
        eqn = exp(sigma_x.*sigma_x./2).*normcdf(sigma_x-(log(x)-sigma_x))-x.*normcdf(-log(x)./sigma_x) == varpi_temp;
        temp = vpasolve(eqn, x);
        if temp>=0
            underline_x_ub(k, i) = temp;
        else
            no_solution(k, i) = 1;
        end
    end
end
%  %underline_x_ub =   [21.6480,   23.4126,   26.1022,   20.9898,   27.3432,   21.7458;
%    22.0538,   24.5482,   23.4180,   24.6617,   25.0367,   22.4453;
%    22.3934,   23.5518,   22.5340 ,  23.2726,   24.8191,   22.2654;
%    20.7960,   23.9193,   21.9694 ,  23.7405 ,  22.1600 ,  20.7735;
%    24.1723 ,  26.8913 ,  23.3186 ,  28.1148 ,  27.5073,  19.4319;
%    22.5726  , 23.4773  , 20.7196  , 23.7081 ,  23.5796,   17.8373]
%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 5 guess big_L x_bar
%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializations for x_bar, big_L
if (FUNCTION_ITERATION == 1) 
    x_bar_Guess = underline_x_ub ./ 10;
    for k = 1:K
        for i = 1:N
            big_L_Guess(k, i) = mu(k, i)*big_L_bar(i);
        end
    end
end
big_L_Guess

x_bar0 = min( max( x_bar_Guess , 0.000001 ) , underline_x_ub - 0.000001 ) ;
big_L0 = big_L_Guess;
x_bar0

% Initialize endogenous variables
x_bar   = x_bar0 ;
big_L   = big_L0 ;

dist_out = 1 ;
iter_out = 1 ;
dist_out_best = 1000000;
while ( dist_out > toler_out && iter_out <= maxiter_out )
    
    % If the maximum number of iterations is reached, we
    % compute all outcomes using the best guess of x_bar / big_L
    if (iter_out == maxiter_out)
        big_L = big_L_best;
        x_bar = x_bar_best;
    end
        
    %%%%%%%
    %STEP 6
    %%%%%%%
    
    % Compute Int_x_bar
    for k = 1:K
        for i = 1:N
            Int_x_bar(k,i) = exp(sigma_x*sigma_x/2)...
                *normcdf(sigma_x-(log(x_bar(k,i))-sigma_x))-x_bar(k,i)*normcdf(-log(x_bar(k,i))/sigma_x);
        end
    end
    Int_x_bar
    
    % Return the value of lognormal CDF at x with 0 mean and sigma variance
    for k = 1:K
        for i = 1:N
            G_x_bar(k,i) = normcdf(log(x_bar(k,i))./sigma_x);
        end
    end
    
    % Compute q_theta
    for k = 1:K
        for i = 1:N
            q_theta(k, i) = min(varpi(k,i)/Int_x_bar(k,i),0.99999);
        end
    end
    q_theta
 
    % Compute theta
    for k = 1:K
        for i = 1:N
            theta(k, i) = ( q_theta(k, i)^(-xsi) - 1 )^(1/xsi);
        end
    end
    theta
    
    % Compute u
    for k = 1:K
        for i = 1:N
            u(k,i) = chi(k,i) / ( theta(k,i)*q_theta(k,i)*(1-G_x_bar(k,i)) + chi(k,i) );
        end
    end
    u
    
    % Compute p
    for k = 1:K
        for i = 1:N
            p(k,i) = 1 - theta(k,i)*q_theta(k,i)*(1-G_x_bar(k,i));
        end
    end
    p
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 7: Compute big_L_tilde
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_L_tilde = big_L.*(1-u).*stdint(x_bar, sigma_x);
    big_L_tilde
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 8: Compute w_tilde
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    w_tilde = (gamma.*big_Y)./ big_L_tilde;
    w_tilde
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 9: Compute big_E_V
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_E_V = kappa_tilde.*w_tilde.*theta.*u.*big_L;
    big_E_V
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 10: Compute big_E_C
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:N
        temp1 = 0;
        temp2 = 0;
        for k = 1:K
            temp1 = temp1 + gamma(k, i)*big_Y(k, i);
            temp2 = temp2 + big_E_V(k,i);
        end
        big_E_C(i) = temp1 - temp2 - NX_Data(i,1);
    end
    big_E_C
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 11: Compute lambda_tilde
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lambda_tilde = big_L_bar./big_E_C;
    lambda_tilde
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 12: Solve Bellman Equation
    % and obtain big_U
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_U = compute_U(0.1, 100);
    if (sum(sum(isnan(big_U))) ~= 0 || sum(sum(isinf(big_U))) ~= 0)
        disp("big_U contains NaN or Inf elements");
        quit;
    end
    big_U
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Step 11: Update big_L
    %%%%%%%%%%%%%%%%%%%%%%
    
    %STEP 11a
    denominator = zeros(K, N);
    for k = 1:K
        for i = 1:N
            for l = 1:K
                denominator(k, i) = denominator(k, i)+exp((-C(k,l,i)+b(i)+theta(l,i)...
                    *kappa_tilde(l,i)*lambda_tilde(i)*w_tilde(l,i)*(beta/(1-beta))...
                    +delta*big_U(l,i))/zeta);
            end
        end
    end
    s = zeros(K,K,N);
    for k = 1:K
        for i = 1:N
            for l = 1:K
                s(k, l, i) = (exp((-C(k,l,i)+b(i)+theta(l,i)...
                    *kappa_tilde(l,i)*lambda_tilde(l,i)*w_tilde(l,i)*(beta/(1-beta))...
                    +delta*big_U(l,i))/zeta))/denominator(k, i);
            end
        end
    end
    
    %STEP 11b
    %   return the null space of s_i.
    I = eye(K);
    for i = 1:N
        s_i = s(:,:,i);
        A = I-s_i';
        y_tilde = zeros(K, N);
        big_L1 = zeros(K, N);
        spower_old = eye(K);
        spower_new = zeros(K, K);
        dist_spower = 1;
        iter_spower = 1;
        while (dist_spower > 0.000001 & iter_spower <= 2000)
            spower_new = A*A;
            dist_spower = spower_new - spower_old;
            spower_old = spower_new;
            iter_spower = iter_spower + 1;
        end
        for k = 1:K
            for n = 1:N
                y_tilde(k, n) = spower_new(k, n)/u(k,n);
            end
        end
        phi0   = big_L_bar(i,1) / sum( y_tilde(:,i) );
        for k = 1:K
            big_L1(k, i) = max( phi0*y_tilde(k, i), 0 );
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Step 14: Update x_bar
    %%%%%%%%%%%%%%%%%%%%%%
    x_bar1 = zeros(K,N);
    for k = 1:K
        for i = 1:N
            x_bar1 = ( (1-delta)*big_U(k,i) - eta(k, i) ) / (lambda_tilde(i)*w_tilde(k,i));
            x_bar1 = max( min( x_bar1 , underline_x_ub-0.000001 ) , 0.000001 );
        end
    end
   
    % Initialize steps
    step_L = step_L_ini;
    step_x = step_x_ini;
   
    
    % Compute the distance for loop control
    dist_L = max( max( abs(( big_L1 - big_L)/(max(abs(big_L),0.0000001)) ) ) );
    dist_x = max( max( abs((x_bar1 - x_bar) /(max(abs(x_bar),0.0000001)) ) ) );
    dist_out = max(dist_x,dist_L);
    dist_out_vec(iter_out) = dist_out;
   
    % Saves x_bar, big_L, x_bar1, and big_L1 
    % leading to the best dist_out
    if (dist_out < dist_out_best)
        dist_out_best = dist_out ;
        x_bar_best    = x_bar ;
        big_L_best    = big_L ;
    end
    flag_osc = 0;
    % if iterations don't decrease dist_out, change the value of step.
    if(iter_out > 300)
        % If dist_out does not improve, or if it improves very little (by less than toler_out/10),
        % reduce step size and go back to x_bar and L that led to 
        % the smallest value of dist_out
        if ( flag_osc == 0 || iter_out >= flag_osc + 300 ) 
            mean_dist1 = sum(dist_out_vec(iter_out-49:iter_out,1))/50; 
            mean_dist2 = sum(dist_out_vec(iter_out-99:iter_out-51,1))/50;  
        end 
        if ( ( (mean_dist1 >= mean_dist2) || ...
             ( (mean_dist1 - mean_dist2 < 0) && abs(mean_dist1 - mean_dist2) < toler_out ) ) && ...
             flag_osc == 0) 
            step_L   = 0.9*step_L;
            step_x   = 0.9*step_x;
            flag_osc = iter_out;
        end
        % But we only allow step sizes to be reduced every 300 iterations
        if ( ( (mean_dist1 >= mean_dist2) || ...
             ( (mean_dist1 - mean_dist2 < 0) && abs(mean_dist1 - mean_dist2) < toler_out ) ) && ...
             iter_out >= flag_osc+300 ) 
            step_L   = 0.9*step_L; 
            step_x   = 0.9*step_x;
            flag_osc = iter_out;
        end
    end
   
    big_L_Temp = (1-step_L)*big_L + step_L*big_L1;
    big_L      = min(max(big_L_Temp,.9*big_L),1.1*big_L);
    for i = 1: N
        big_L(:,i) = big_L(:,i)/sum(big_L(:, i))*big_L_bar(i,1);
    end
    x_bar = (1-step_x)*x_bar + step_x*x_bar1 ;
    x_bar = max( min( x_bar , underline_x_ub - 0.000001 ) , 0.000001 );
    
    % If step size was adjusted after bad behavior, 
    % restart with x_bar = x_bar_best 
    % and big_L = big_L_best instead   
    if (flag_osc == iter_out)
        big_L = big_L_best ;
        x_bar = x_bar_best ;
    end
    
    iter_out = iter_out + 1;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Simulated Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Emp_Count = zeros(K, N);
ShareEmp = zeros(K, N);
Unemp = zeros(N);

%%%%%%%%%%%%%%%%%%%
% Employment Shares
%%%%%%%%%%%%%%%%%%%

for k = 1:K
    for i = 1:K
        Emp_Count(k,i) = big_L(k, i)*(1-u(k,i)) ;
        ShareEmp(k, i)  = Emp_Count(k,i) / sum(Emp_Count(:, i)) ;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average unemployment by country
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
    term1 = 0;
    term2 = 0;
    for k = 1:K
        term1 = term1 + big_L(k, i)*u(k, i);
        term2 = term2 + big_L(k, i);
    end
    Unemp(i) = term1/term2;
end


%%%%%%%%%%%%%%%
% Average Wages
%%%%%%%%%%%%%%%

for k = 1:K
    for i = 1:N
        term1 = normcdf(sigma_x - log(x_bar(k,i))/sigma_x);
        term2 = normcdf(-log(x_bar(k,i))/sigma_x);
        w_bar(k,i) = (1-beta(k,i))*w_tilde(k,i)*x_bar(k,i) + beta(k,i)*w_tilde(k,i)* ...
        exp( (sigma_x^2)/2 )*term1/term2;
    end
end


%%%%%%%%%%%%%%%%%%%
% Variance of Wages
%%%%%%%%%%%%%%%%%%%

for k = 1:K
    for i = 1:N
        term1 = normcdf(sigma_x - log(x_bar(k,i))/sigma_x);
        term2 = normcdf(-log(x_bar(k,i))/sigma_x);
        term3 = normcdf(2*sigma_x - log(x_bar(k,i))/sigma_x);
        Var_w(k,i) = ( ( beta(k,i)*w_tilde(k,i) )^2 )* ...
        ( exp(2*sigma_x^2)*term3/term2 - ...
          exp(sigma_x^2)*(term1/term2)^2 );
    end
end



%%%%%%%%%%%%%
% Labor Share
%%%%%%%%%%%%%
for k = 1:K
    for i = 1:N
        labor_share(k, i) = ( w_bar(k, i)*Emp_Count(k, i) )/ big_Y(k,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-period transition rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%



s_tilde = zeros(K,K,N);

for i = 1:N
    % unemployment to unemployment
    for l = 1:K
        for k = 1:K
            s_tilde(l,k,i) = s(l,k,i)*(1 - theta(k,i)*q_theta(k,i)*(1-G_x_bar(k,i)));
        end
    end
    
    % unemployment to employment
    for l = 1:K
        for k = 1:K
            s_tilde(l,K+k,i) = s(l,k,i)*theta(k,i)*q_theta(k,i)*(1-G_x_bar(k,i));
        end
    end
    
    % employment to employment
    for l = 1:K
        for k = 1:K
            if ( l == k )
                s_tilde(K+l,K+k, i) = (1-chi(l, i));
            else
                s_tilde(K+l,K+k, i) = 0;
            end
        end
    end
    
    % employment to unemployment
    for l = 1:K
        for k = 1:K
            if ( l == k )
                s_tilde(K+l,k, i) = chi(l,i);
            else
                s_tilde(K+l,k,i) = 0;
            end
        end
    end
end
FUNCTION_ITERATION = FUNCTION_ITERATION+1;
FUNCTION_ITERATION
end
    





function [U] = compute_U(init_lb, init_ub)
%STEP 12

U_old = initiate_U(init_lb, init_ub);
U_new = update_U(U_old);
while(distance(U_old, U_new) > 0.001)
    U_old = U_new;
    U_new = update_U(U_old);
end
U = U_new;
end

function [U] = initiate_U(init_lb, init_ub)
%STEP 12a

global K N;
gap = init_ub -  init_lb;
U = rand(K, N)*gap+init_lb;
end

function [U_new] = update_U(U_old)
%STEP 12b

global K N;
global C b theta kappa_tilde lambda_tilde w_tilde beta delta zeta;
bracket = zeros(K, N) ;
for k = 1:K
    for i = 1:N
        for l = 1:K
            bracket(k, i) = bracket(k, i)+exp((-C(k,l,i)+b(i)+theta(l,i)...
                *kappa_tilde(l,i)*lambda_tilde(i)*w_tilde(l,i)*(beta/(1-beta))...
                +delta*U_old(l,i)-delta*U_old(k,i))/zeta);
        end
    end
end
U_new = zeta.*log(bracket) + delta*U_old;
end

function [distance] = distance(U_old, U_new)
%return the average L2 distance of each elements between U_old and U_new
global K N;
dist = U_old - U_new;
dist_sqrt = dist.*dist;
distance = sum(sum(dist_sqrt))/(N*K);
end


    
    
    
    
    
    
    
   