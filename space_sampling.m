
clear all
clc


%% Sampling the state space

m=input('Enter the number of prices:')
N=input('Enter the parking capacity:')

%m=3;
%N=6;

%%

load D_N50_m4.mat

%% Sample Generation
T=30;
Time_step=1;
n32=1;
sample_length=100;

while Time_step<=T
while n32<= sample_length
    g=(N-1)*rand(m,1);
    g=round(g);
    [ok] = state_difference(g',D(1:m,Time_step));
    if sum(g)<=N & ok==1
        sample_space(n32,1:m,Time_step)=g';
        n32=n32+1;
    end
end
Time_step=Time_step+1;
n32=1;
end

%%
%r_hat =[725;3.84;5.02;2.55];
load r_hat_N50_m4.mat

lambda_max=[0.55;0.5;0.3;0.2];
mu_min=[0.2;0.2;0.4;0.45];
lambda_min=[0.3;0.2;0.1;0.08];
mu_max=[0.5;0.6;0.6;0.62];

%lambda=[0.6;0.5;0.3;0.2];
%mu=[0.2;0.2;0.4;0.4];
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
    k=k-1
    %feauters=[constant sample_sapce(:,:,k)];
end


