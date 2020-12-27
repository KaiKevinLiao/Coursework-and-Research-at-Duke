function y = LiFunc_3 (para,numP,numI,L,Z,X,I) 
    
    alpha = para(1);
    gamma = para(2:end);
    
    theta_A = zeros(numP,1);
    theta = ones(numP,1);
    
    sigma = X(6,:);            
    s = zeros(1,numP);          
    expm = zeros(numI,numP);    
    denom = zeros(numI,1);     
    ind_llh = zeros(numI,numP); 
    
    %%%%%%%%%%%%%
    % Contraction Mapping
    %%%%%%%%%%%%%
    
    while sum((theta - theta_A).^2) > 0.00001  %while theta not converge
        
        theta = theta_A; % update theta
        for i = 1:numI         
            for j = 1:numP
                expm(i,j) = exp(theta(j,1) + alpha * Z(i,j) +...,
                    gamma(1) * X(5,j) * I(2,i)+..., %panfish with kids
                    gamma(2) * X(2,j) * I(2,i)+..., %restroom with kids
                    gamma(3) * X(1,j) * I(1,i)+..., %ramp with boat
                    gamma(4) * X(3,j) * I(1,i));    %walleye with boat
            end
        end
        for i = 1:numI
            denom(i,1) = sum(expm(i,:));
            for j = 1:numP
                ind_llh(i,j) = expm(i,j) / denom(i,1); %the likelihood of person i going to place j
            end
        end
        s = sum(ind_llh) / numI;        % predicted share for each place
        for j = 1:numP
            theta_A(j,1) = theta(j,1) + (log(sigma(1,j)) - log(s(1,j)));    %update theta one time(Berry Contraction mapping)       
        end
    end
    mean_theta = mean(theta_A);
    theta_A = theta_A - mean_theta;  % demean
    
    %%%%%%%%%%%%%
    % Likelihood Function
    %%%%%%%%%%%%%

	LM = zeros(numI,1); %likelihood matrix: the i-th row is the probablity of individual i choose j and not others 
    expm = zeros(numI,100);     %exponential matrix
	for i = 1:numI
        c = L(i,1);  % choice s
        for j = 1:numP
            expm(i,j) = exp(theta_A(j,1) + alpha * Z(i,j)+...,
                gamma(1) * X(5,j) * I(2,i)+..., %panfish with kids
                gamma(2) * X(2,j) * I(2,i)+..., %restroom with kids
                gamma(3) * X(1,j) * I(1,i)+..., %ramp with boat
                gamma(4) * X(3,j) * I(1,i));    %walleye with boat); %the "raw" likelihood of person i going to place j
        end
     
        LM(i,1) = log(expm(i,c)) - log(sum(expm(i,:))) ; %likelihood of a single person i
	end
	y = -sum(LM);
end