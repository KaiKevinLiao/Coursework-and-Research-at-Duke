function y = LiFunc_2 (para,numP,numI,L,Z,X,I) 

    alpha = para(1);
    beta = para(2:6);
    gamma = para(7:10);

    % likelihood for each individual choosing the correct site
	likelihood = zeros(numI,1);  
    expm = zeros(numI,numP); 
	for i = 1:numI
        c = L(i,1);  % choice 
        for j = 1:numP
            expm(i,j) = exp(beta(1) * X(1,j) + beta(2) * X(2,j) + beta(3) * X(3,j) ...
                + beta(4) *X(4,j) + beta(5) *X(5,j) + gamma(1) * X(5,j) * I(2,i) +...,
                gamma(2) * X(2,j) * I(2,i) + gamma(3) * X(1,j) * I(1,i) + gamma(4) * X(3,j) *I(1,i) ...
                + alpha * Z(i,j)); 
        end     
        likelihood(i,1) = log(expm(i,c)) - log(sum(expm(i,:))) ; 
	end
	y = -sum(likelihood);
end
