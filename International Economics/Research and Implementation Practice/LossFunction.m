function [F] = LossFunction(param_in)
%LossFunction Compute the loss function of estimation

global gamma nu mu pi NX_Data big_Y_Data Unemp_Data Emp_Data ...
    ShareEmp_Data big_L_bar_Data w_bar_Data CoeffVar_w_Data s_tilde_Data...
    Ergodic_Dist_Data big_L_bar World_big_Y_Data NX_Data;
global COUNTRY_EST K N numP NP NP_fixed TempSimAnn FUNCTION_ITERATION;
global s s_tilde Unemp ShareEmp w_bar Var_w  labor_share big_L  lambda_tilde;
global param_in;
global x_bar w_tilde q_theta theta G_x_bar dist_out;
global b eta kappa_tilde sigma_x chi C_US C;

global lambda delta xsi zeta beta;

global big_Y;

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

ComputeEQ()

WeightU = 1;
WeightSh = 1;
WeightW = 1;
WeightVacCosts = 1;

F_U = 0;
F_Sh = 0;
F_W = 0;
F_VacCosts = 0;



VacancyCosts_Data = 0.58; % why?
% compute VacancyCosts from EQ
for k = 1:K
    for i = 1:N
        VacancyCosts(k,i) = kappa_tilde(k, i)*w_tilde(k, i) / w_bar(k, i);
    end
end

for i = 1:N
    F_U = F_U + WeightU*abs( ( Unemp(i) - Unemp_Data(i) )...
        / (abs(Unemp_Data(i))) );
end

for i = 1:N
    for k = 1:K
        F_sh = F_sh + WeightSh*abs( (ShareEmp(k,i) - ShareEmp_Data(k,i))...
            / (abs(ShareEmp_Data(k,i))) );
    end
end

for i = 1:N
    for k = 1:K
        F_W = F_W + WeightW*abs( (w_bar(k,i) - w_bar_Data(k,i))...
            / (abs(w_bar_Data(k,i)) ) );
    end
end

for i = 1:N
    for k = 1:K
        F_VacCosts = F_VacCosts + WeightVacCosts*...
            abs( WeightVacCosts*abs( VacancyCosts(k, i) - VacancyCosts_Data )...
            / abs(VacancyCosts_Data) );
    end
end

toler_out = 0.0005;
F_p1 = 1000000.0*max(dist_out - toler_out,0.0) ;


F = F_Tr + F_U + F_Sh + F_W  + F_VacCosts + ...
        F_p1;

end

