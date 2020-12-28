function [x_new] = update_x_underline(x, w_tilde, U, eta, delta, lambda_tilde, N, K)
%STEP14
%   This question is not finished due to a question
lambda_expand = repmat(lambda_tilde,K,1);
x_pre = zeros(K,N);
x_pre = ((1-delta)*U-eta)./(lambda_expand.*w_tilde);
x_new = min(x_pre, )
end

