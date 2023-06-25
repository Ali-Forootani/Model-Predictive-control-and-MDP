function [lambda mu] = lambda_mu_calculation(s,N,m,lambda_max,lambda_min,mu_max,mu_min)

for i=1:m
    %mu(i,1)=mu_max(i)+(mu_min(i)-mu_max(i))*exp(-s(i)/N);
    mu(i,1)=mu_max(i)+(mu_min(i)-mu_max(i))*exp(-s(i)/N);
    
    %mu(i,1)=mu_max(i)*(1-exp(-s(i)/N));
    
    lambda(i,1)=lambda_min(i)+(lambda_max(i)-lambda_min(i))*exp(-s(i)/N);
end

