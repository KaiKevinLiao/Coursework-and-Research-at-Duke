function y = LiFunc_1 (para,numP,numI,L,Z,X) %likelihood function

    alpha = para(1);
    beta = para(2:end);

    % likelihood for each individual choosing the correct site
	likelihood = zeros(numI,1);  
    expm = zeros(numI,numP); 
	for i = 1:numI 
        for j = 1:numP
            expm(i,j) = exp(beta(1) * X(1,j) + beta(2) * X(2,j) + beta(3) * X(3,j) ...
                + beta(4) *X(4,j) + beta(5) *X(5,j) + beta(6) * X(6,j) + alpha * Z(i,j)); 
        end     
        c = L(i,1); 
        likelihood(i,1) = log(expm(i,c)) - log(sum(expm(i,:))) ; 
	end
	y = -sum(likelihood);
end