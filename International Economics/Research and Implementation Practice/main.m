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

FUNCTION_ITERATION = 0;
set_global_data();
read_data();
read_parameters();
%Compute B matrix
B = zeros(N*K, N*K);
for k = 1:K
    for o = 1:N
        for l = 1:K
            for i = 1:N
                B((k-1)*N+o, (l-1)*N+i) = pi(k, o,i)*(mu(k,i)*gamma(l,i)+...
                    (1-gamma(l,i))*nu(l,k,i));
            end
        end
    end
end

%Compute A matrix
A = [eye(N*K) - B; ones(1, N*K)];

%Compute b matrix
b_rhs = zeros(N*K+1, 1);
NX_Data(1)
for k = 1:K
    for o = 1:N
        temp = 0;
        for i = 1:N
            temp = temp + pi(k, o, i)*mu(k, i)*NX_Data(i);
        end
        b_rhs((k-1)*N+o,1) = -temp;
    end
end
b_rhs(N*K+1, 1) = 1;

%Solve for Y by AY = b
Y = linsolve(A, b_rhs);

%Write big_Y
for k = 1:K
    for o = 1:N
        big_Y(k, o) = Y((k-1)*N+o);
    end
end

ObjectiveFunction = @LossFunction;
rng default % For reproducibility
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,param_in);






