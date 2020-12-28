function [omegabar] = get_omegabar(chi, delta, kappatilde, beta)
%STEP3.1 get the value of omegabar
%   
omegabar = ((1-(1-chi)*delta).*kappatilde)./(delta*(1-beta));
end

