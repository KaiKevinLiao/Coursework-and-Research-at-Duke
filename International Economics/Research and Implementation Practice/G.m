function [g] = get_G(x,sigma)
%LOGNORMAL CDF 
%   Return the value of lognormal CDF at x with 0 mean and sigma variance
g = normcdf(log(x)./sigma);
end

