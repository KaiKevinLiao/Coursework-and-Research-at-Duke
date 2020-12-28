function [Int] = stdint(x, sigma)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Int = exp(sigma.*sigma./2).*(normcdf(sigma-(log(x)./sigma))./normcdf(-(log(x)./sigma)));
end

