function [I] = compute_I(x,sigma)
%STEP 6
if x == 0
    I = exp(sigma.*sigma./2);
else
    I = exp(sigma.*sigma./2).*normcdf(sigma-(log(x)-sigma))-x.*normcdf(-log(x)./sigma);
end
end

function [g] = compute_G(x, sigma) 
%STEP 6
%   Return the value of lognormal CDF at x with 0 mean and sigma variance
g = normcdf(log(x)./sigma);
end

function [theta] = compute_theta(omega_bar, I, xi, N, K)
%STEP 6
%   we assume xi is an array from 1 to N
xi_expand = repmat(xi,K,1); %expand xi to an K*N matrix
y = omega_bar./I;
theta = ((1-y.^xi)./(y.^xi)).^(1.^xi);
end

function [u] = compute_u(chi, theta, x, sigma, xi, N, K) 
%STEP 6
xi_expand = repmat(xi,K,1); %expand xi to an K*N matrix
q = ((1-theta.^xi)./(theta.^xi)).^(1.^xi);
g = compute_G(x, sigma);
u = chi./(theta.*q.*(1-g))+chi;
end

function[L_tilde] = compute_L_tilde(L, u, x, sigma)
%STEP 7
L = L.*(1-u).*stdint(x, sigma);
end

function[w_tilde] = compute_L_tilde(gamma, Y, L_tilde)
%STEP 8
w_tilde = (gamma.*Y)./L_tilde;
end

function[E_v] = compute_E_v(kappa_tilde, w_tilde, theta, u, L)
%STEP 9
E_v = kappa_tilde.*w_tilde.*theta.*u.*L;
end

function[E_c] = compute_E_c(gamma, Y, E_v, NX)
%STEP 10
E_c = sum(gamma.*Y) - sum(E_v) - NX;
end

function[lambda_tilde] = compute_lambda_tilde(L, E_c)
%STEP 11
lambda_tilde = sum(L)./E_c;
end