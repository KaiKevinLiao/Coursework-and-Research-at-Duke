function [I] = get_I(x,sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if x == 0
    I = exp(sigma.*sigma./2);
else
    I = exp(sigma.*sigma./2).*normcdf(sigma-(log(x)-sigma))-x.*normcdf(-log(x)./sigma);
end
end

