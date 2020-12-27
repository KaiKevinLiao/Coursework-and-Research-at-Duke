function y = LiFunc_6 (para,numP,numI,numK,L,Z,X,I) 
    
    alpha = para(1);
    gamma = para(2:5);
    sigma_v = para(6);

    theta_A = zeros(numP,1);
    theta = ones(numP,1);
    
    sigma = X(6,:);             % Data Share of Each Site
    s = zeros(1,numP);          % Predicted Share of Each Site
    expm = zeros(numI,numP);    % Matrix for Exp
    ind_llh = zeros(numI,numP); % Likelihood Matrix
    
    exp_temp = zeros(numK,numP);
    llh_temp = zeros(numK,numP);
    denom_temp = zeros(numK,1);
    v = randn(numI,numK) * sigma_v;    % random beta
    
    %%%%%%%%%%
    % Compute Theta by Contraction Mapping
    %%%%%%%%%%
    
    while sum((theta - theta_A).^2) > 0.00001  
        
        theta = theta_A; % update theta
        for i = 1:numI  
            for k = 1:numK
                for j = 1:numP
                    exp_temp(k,j) = exp(theta(j,1) + alpha * Z(i,j) +...,
                        gamma(1) * X(5,j) * I(2,i)+..., 
                        gamma(2) * X(2,j) * I(2,i)+..., 
                        gamma(3) * X(1,j) * I(1,i)+..., 
                        gamma(4) * X(3,j) * I(1,i)+..., 
                        v(i,k)   *  X(3,j));    
                end
            end
            denom = sum(exp_temp,2);
            for k = 1:numK      
                for j = 1:numP
                    llh_temp(k,j) = exp_temp(k,j) / denom(k,1);
                end
            end
            ind_llh(i,:) = mean(llh_temp);
        end
        
        s = sum(ind_llh) / numI;        
        
        % Update theta
        
        for j = 1:numP
            theta_A(j,1) = theta(j,1) + (log(sigma(1,j)) - log(s(1,j)));           
        end
        
    end
    % Normalization of Theta
    mean_theta = mean(theta_A);
    theta_A = theta_A - mean_theta; 
    
    %%%%%%%%%%%%%%%%
    % Compute Likelihood Using Theta_A
    %%%%%%%%%%%%%%%%
    
    LM_temp = zeros(numK, numP);
    LM = zeros(numI,1);
    for i = 1:numI 
        c = L(i,1);
        for k = 1:numK
            for j = 1:numP
                exp_temp(k,j) = exp(theta_A(j,1) + alpha * Z(i,j) +...,
                    gamma(1) * X(5,j) * I(2,i)+..., 
                    gamma(2) * X(2,j) * I(2,i)+..., 
                    gamma(3) * X(1,j) * I(1,i)+..., 
                    gamma(4) * X(3,j) * I(1,i)+..., 
                    v(i,k)   *  X(3,j));    
            end
        end
        denom = sum(exp_temp,2);
        for k = 1:numK      
            for j = 1:numP
                LM_temp(k,j) = log(exp_temp(k,j)) - log(denom(k,1));
            end
        end
        LM(i,1) = mean(LM_temp(:,c));
    end
    y = -sum(LM);
    
end