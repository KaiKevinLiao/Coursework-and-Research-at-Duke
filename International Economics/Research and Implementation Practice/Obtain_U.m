function [U] = compute_U(zeta, C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, delta, N, K, threshold, init_lb, init_ub)
%STEP 12
%   Inputs include a threshold of convergence. It is a float number. 
%   The first guess U is drawn from a uniformly distribution between
%   init_lb and init_ub.
U_old = initiate_U(N, K, init_lb, init_ub);
U_new = update_U(U_old, zeta, C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, delta, N, K);
while(distance(U_old, U_new, N, K) > threshold)
    U_old = U_new;
    U_new = update_U(U_old, zeta, C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, delta, N, K);
end
U = U_new;
end

function [U] = initiate_U(N, K, init_lb, init_ub)
%STEP 12a
gap = init_lb, init_ub;
U = rand(K, N)*gap+init_lb;
end

function [U_new] = (U_old, zeta, C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, delta, N, K)
%STEP 12b
%   need to check the format of C and zeta. Here we assume zeta is an array
%   and C is a k*l*i matrix
%   This code does not use many matrix operation. May be slow.
zeta_expand = repmat(zeta,K,1); %expand xi to an K*N matrix
%Calculate the elements in the blanket
bracket = zeros(K, N) ;
for k = 1:K
    for i = 1:N
        for l = 1:K
            bracket(k, i) = bracket(k, i)+exp((-C(k,l,i)+b(l,i)+theta(l,i)...
                *kappa_tilde(l,i)*lambda_tilde(l,i)*w_tilde(l,i)*(beta(l,i)/(1-beta(l,i)))...
                +delta*U_old(l,i)-delta*U_old(k,i))/zeta_expand(k, i));
        end
    end
end
U_new = zeta_expand.*log(bracket) + delta*U_old;
end

function [distance] = distance(U_old, U_new, N, K)
%return the average L2 distance of each elements between U_old and U_new
dist = U_old - U_new;
dist_sqrt = dist.*dist;
distance = sum(sum(dist_sqrt))/(N*K);
end

