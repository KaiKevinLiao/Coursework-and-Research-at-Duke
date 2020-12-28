function read_parameters

%%%%%%%%%%
% Read Starting_Point.csv
%%%%%%%%%%

cd 'C:\Users\kaike\OneDrive\Desktop\International Trade\Dix-Carneiro 2020 Code'
global K N param_in
global lambda delta xsi zeta beta param_lb param_ub;

% Read Table
temp_storage = readtable('Starting_Point.csv');
temp = temp_storage(:, 1:6);
temp = table2array(temp);
param_in = temp(:,2);
param_lb = temp(:,4);
param_ub = temp(:,5);


global b eta kappa_tilde sigma_x chi C_US C;

% Read b
for n = 1:N
    b(n) = temp(n, 2);
end

% Read kappa_tilde
for k = 1:K
    for n = 1:N
        kappa_tilde(k, n) = temp(6+(n-1)*K+k, 2);
    end
end

% Read sigma_x
sigma_x = temp(43, 2);

% Read chi
% chi(sector, country) = chi(country) + chi(sector)
chi_country = zeros(N);
for n = 1:N
    chi_country(n) = temp(43+n, 2);
end
% Read chi_sector
chi_sector = zeros(K);
for k = 1:K
    chi_sector(k) = temp(49+n, 2);
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
        eta(k, n) = temp(55+(n-1)*K+k, 2);
    end
end

% Read C_US
for k = 1:K
    for l = 1:K
        C_US(k,l) = temp(91+(k-1)*K+l, 2);
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

clear temp_storage temp;

temp_storage = readtable('FixedParams.csv');
temp = temp_storage(:, 1:2);
temp = table2array(temp);
FixedPar = temp(:,2);
lambda = FixedPar(1);
delta  = FixedPar(2); 
xsi    = FixedPar(3);
zeta   = FixedPar(4);
beta   = FixedPar(5);


end

