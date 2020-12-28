function [L] = update_L(C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, U, delta, zeta, u, N, K)
%STEP 13
s = compute_s(C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, U, delta, zeta, N, K);
y = compute_y(s, N, K);
L = compute_L(u, y);
end

function [s] = compute_s(C, b, theta, kappa_tilde, lambda_tilde, w_tilde, beta, U, delta, zeta, N, K)
%STEP 13a
denominator = zero(K, N);
for k = 1:K
    for i = 1:N
        for l = 1:K
        denominator(k, i) = denominator(k, i)+exp((-C(k,l,i)+b(l,i)+theta(l,i)...
                *kappa_tilde(l,i)*lambda_tilde(l,i)*w_tilde(l,i)*(beta(l,i)/(1-beta(l,i)))...
                +delta*U(l,i))/zeta(i));
        end
    end
end
s = zeros(K,K,N);
for k = 1:K
    for i = 1:N
        for l = 1:K
            s(k, l, i) = (exp((-C(k,l,i)+b(l,i)+theta(l,i)...
                *kappa_tilde(l,i)*lambda_tilde(l,i)*w_tilde(l,i)*(beta(l,i)/(1-beta(l,i)))...
                +delta*U(l,i))/zeta(i)))/denominator(k, i);
        end
    end
end
end

function [y] = compute_y(s, N, K)
%STEP 13b
%   return the null space of s_i.
I = eye(K);
for i = 1:N
    s_i = s(:,:,i);
    A = I-s_i';
    y(:,i) = null(s_i);
end
end

function [L] = compute_L(u, y)
%STEP 13c
end