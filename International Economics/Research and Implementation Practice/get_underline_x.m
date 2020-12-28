function [underline_x_ub, no_solution] = get_underline_x (omega_bar, sigma, N, K)
%STEP 4 
%   return the upperbound for x, if not, return the position where this no
%   feasiable numerical solution
no_solution = zeros(K,N);
for k = 1:K
    for i = 1:N
        syms x;
        sigma = sigma(k, i);
        omega = omega_bar(k, i);
        eqn = exp(sigma.*sigma./2).*normcdf(sigma-(log(x)-sigma))-x.*normcdf(-log(x)./sigma) == omega;
        s = vpasolve(eqn, x);
        if x>=0
            underline_x_ub(k, i) = s;
        else
            no_solution(k, i) = 1;
    end
end
end

 