function [t,a,b] = check_omega_bar(omega_bar, sigma, N, K)
%STEP 3 Check if omega_bar satisfies q(theta) <= 0
%   Detailed explanation goes here
a = [];
b = [];
t = [];
for k = 1:K
    for i = 1:N
        if omega_bar(k, i)/stdint(0,sigma(k, i)) > 1
            a = [a, k];
            b = [b, i];
            t = [t, 1];
            % this is yet to finished, we need to include how to abort the
            % procedure or penalize the objective function
        else
            t = [t, 0];
        end
    end
end
end

