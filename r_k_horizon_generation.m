function [r_k_horizon]= r_k_horizon_generation(m,N,sample_space,sample_length,lambda_max,lambda_min,mu_max,mu_min,C,D,T)

load r_hat_N50_m4.mat
lambda=[0.6;0.5;0.3;0.2];
mu=[0.2;0.2;0.4;0.4];
C=[0.9 1 1.1 1.2];
constant=ones(sample_length,1);
%%

k=T-1;
d=T;
r_k(1:m+1,d)=r_hat;
%n1=
while k >= 1
    [Revenue] = DP_Space_sampling(m,N,sample_space(:,:,k),sample_space(:,:,k+1)...
    ,D,lambda_max,lambda_min,mu_max,mu_min,T,k+1,C,r_hat);
    features=[constant sample_space(:,:,k)];
    values=Revenue(:,m+1,1);
    %r_k(:,k)=feauters(:,1:m+1)\values;
    %r_k(1:m+1,d)=r_hat;
    beta=1500*eye(m+1);
    f=features'*features+beta;
    f_1=inv(f);
    gamma=features'*values+1500*r_hat;
    r_hat=f_1*gamma;
    %r_hat=f_2*values;
    %r_hat=features(:,1:m+1)\values;
    r_k(1:m+1,k)=r_hat;
    d=d-1;
    k=k-1;
    %feauters=[constant sample_sapce(:,:,k)];
end
%%
r_k_horizon=r_k;