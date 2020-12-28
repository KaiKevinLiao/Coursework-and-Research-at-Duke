function read_data

%%%%%%%%%%
% Read data
%%%%%%%%%%

cd 'C:\Users\kaike\OneDrive\Desktop\International Trade\Dix-Carneiro 2020 Code'
global K N;
global  gamma nu mu pi NX_Data big_Y_Data Unemp_Data Emp_Data ...
    ShareEmp_Data big_L_bar_Data w_bar_Data CoeffVar_w_Data s_tilde_Data ...
    Ergodic_Dist_Data big_L_bar World_big_Y_Data NX_Data

% Read gamma
temp_storage = readtable('gammaMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global gamma;
gamma = zeros(K, N);
for n = 1:N
    for k = 1:K
        gamma(k, n) = temp((n-1)*K+k, 4);
    end
end
clear temp_storage temp;

% Read nu
% The format is nu(User, Supplier, Country)
temp_storage = readtable('nuMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global nu;
nu = zeros(K, K, N);
for n = 1:N
    for k = 1:K
        for l = 1:K 
            nu(k, l, n) = temp((n-1)*K*K+(k-1)*K+l, 5);
        end
    end
end
clear temp_storage temp;

% Read mu
temp_storage = readtable('muMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global mu;
mu = zeros(K, N);
for n = 1:N
    for k = 1:K
        mu(k, n) = temp((n-1)*K+k, 4);
    end
end
clear temp_storage temp;

% Read pi
%The format is nu(Sector, Exporter, Importer)
temp_storage = readtable('piMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global pi;
pi = zeros(K, N, N);
for k = 1:K
    for n = 1:N
        for l = 1:N 
            pi(k, n, l) = temp((n-1)*N*N+(l-1)*K+k, 5);
        end
    end
end
clear temp_storage temp;

% Read NX_Data (Net Export)
temp_storage = readtable('nxMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global nx;
nx = zeros(N);
for n = 1:N
    nx(n) = temp(n, 3);
end
clear temp_storage temp;

% Read big_Y_Data
temp_storage = readtable('goMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global big_Y_Data;
big_Y_Data = zeros(K, N);
for n = 1:N
    for k = 1:K
        big_Y_Data(k, n) = temp((n-1)*K+k, 4);
    end
end
clear temp_storage temp;

% Read Unemp_Data (Unemployment data)
% Only have data in country level, not secter specific. Measure by
% proportion. 
temp_storage = readtable('unemp_ILO.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global Unemp_Data;
Unemp_Data = zeros(N);
for n = 1:N
    Unemp_Data(n) = temp(n, 3);
end
% rescale
Unemp_Data = Unemp_Data ./ 100;
clear temp_storage temp;

% Read Emp_Data
% Sector specific. Measured by the THOUSANDS of number of employment
% workers. Rescale to the exact numbers.
temp_storage = readtable('goMat.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global Emp_Data;
Emp_Data = zeros(K, N);
for n = 1:N
    for k = 1:K
        Emp_Data(k, n) = temp((n-1)*K+k, 4);
    end
end
% rescale
Emp_Data = Emp_Data * 1000;
clear temp_storage temp;

% Construct ShareEmp_Data
% Share of employment in a specific sector for each country. 
global ShareEmp_Data;
ShareEmp_Data = zeros(K, N);
temp = sum(Emp_Data);
for k = 1:K
    for n = 1:N
        ShareEmp_Data(k, n) = Emp_Data(k, n)/temp(n);
    end
end
ShareEmp_Data
clear temp;

% Construct big_L_bar_Data
% Total number of worker in a country.
global big_L_bar_Data;

temp = sum(Emp_Data);
for n = 1:N
    big_L_bar_Data(n) = temp(n)/(1-Unemp_Data(n));
end
clear temp;

% Read w_bar_Data
% Wage
temp_storage = readtable('wiot_sea_Wage.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global w_bar_Data;
w_bar_Data = zeros(K, N);
for n = 1:N
    for k = 1:K
        w_bar_Data(k, n) = temp((n-1)*K+k, 4);
    end
end
clear temp_storage temp;

% Read CoeffVar_w_Data
% Specific in sectors
temp_storage = readtable('CPS_CoefVar.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global CoeffVar_w_Data;
CoeffVar_w_Data = zeros(K, N);
for k = 1:K
    CoeffVar_w_Data(n) = temp(k, 3);
end
clear temp_storage temp;

% Read s_tilde_Data
% Transition. K+1 by K+1 transition matrix, s_tilde_Data(origin, distination)
temp_storage = readtable('CPS_transitions.csv', 'HeaderLines',1);
temp = table2array(temp_storage);
global s_tilde_Data;
s_tilde_Data = zeros(K+1, K+1);
for n = 1:K+1
    for k = 1:K+1
        s_tilde_Data(n, k) = temp((n-1)*(K+1)+k, 4);
    end
end
A = sum(s_tilde_Data, 2);
% Make sure rows of s_tilde_Data sum exactly to 1
for n = 1:K+1
    for k = 1:K+1
        s_tilde_Data(n, k) = s_tilde_Data(n, k)/A(n);
    end
end
clear temp_storage temp;

% Ergodic Distribution Implied by s_tilde_Data
Ergodic_Matrix_old = eye(K+1);
Ergodic_Matrix_new = zeros(K+1, K+1);
max_dist_erg = 1;
while(max_dist_erg > 0.000001)
    Ergodic_Matrix_new = Ergodic_Matrix_old * s_tilde_Data;
    max_dist_erg = abs(max(max(Ergodic_Matrix_new - Ergodic_Matrix_old)));
    Ergodic_Matrix_old = Ergodic_Matrix_new;
end
% There is a unsolved problem here:
% Ergodic_Dist_Data is very different from previous ShareEmp_Data, it seems
% that it is reversed (sector 1 should be the lowest, not highest). Need to
% check with Fortran output
Ergodic_Matrix_new
global Ergodic_Dist_Data;
Ergodic_Dist_Data = transpose(Ergodic_Matrix_old(1:1,:));
Ergodic_Dist_Data
Unemp_Data(1,1)    = Ergodic_Dist_Data(1,1);
ShareEmp_Data(:,1) = Ergodic_Dist_Data(2:K+1,1) / (1-Unemp_Data(1,1));
ShareEmp_Data

% Normalize country sizes -- US has size 1
L_US = big_L_bar_Data(1);
global big_L_bar;
big_L_bar = big_L_bar_Data / L_US;

% World Gross Output is normalized to 1
% Normalize other nominal variables accordingly
big_Y_Data
global World_big_Y_Data NX_Data;
World_big_Y_Data = sum( sum(big_Y_Data)  );
w_bar_Data = (w_bar_Data / World_big_Y_Data)*L_US;
NX_Data = nx / World_big_Y_Data;
big_Y_Data = big_Y_Data / World_big_Y_Data;
big_Y_Data
sum(sum(big_Y_Data))

end