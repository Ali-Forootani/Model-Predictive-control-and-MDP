function [next_state] = next_state_generation_ADP_Booking_unbounded_bounded_12(state,N,m,D_k,...
D_k_p,F_k,F_k_p,lambda_max,lambda_min,mu_max,mu_min,C,Decision)
%%

f=rand;

less_less=0;
equal_less=0;
curr_demand=0;
full_system=0;





%================ Exact DP Booking Easter Period 

%clear all
%clc
%tic
%%
%m=input('Enter the number of prices:')
%N=input('Enter the parking capacity:')
%m=3;
%N=6;


%% Time Horizon, price vector,
% T: variable for the horizon
%  C: Price vector

%D_k=[0;0;0];
%D_k_p=[1;1;1];
%lambda=[0.6;0.5;0.3];
%mu=[0.2;0.2;0.4];
%lambda=[0.6;0.5];
%mu=[0.2;0.2];
%T=30;
%C=[0.9 1];
%C=[0.9 1 1.1];
%% State Space Construction
% S: state_space
% base(m,N,n): is a written function which get as inputs ,(m,N,a number in
% decimal) and give a possible candid for being a state of the system
% we want to know in how many ways "a decimal number" could be divided
% among "m" different groups!

% n=0;
% l=0;
% %m=3;
% %N=4;
% while n<=((N+1)^m - 1)
%     s=base(m,N,n);
%     if sum(s)<= N
%     l=l+1;
%     S(l,1:m)=s;
%     end
%     n=n+1;
% end

% NS: size of entire state space
%NS=size(S);


%% Generating the reservation vector
% 
%n32=1;
%T=60;
%m=3;
% res_norm=3;
% n=0;
% l=0;
% while n<=((res_norm+1)^m - 1)
%     s=base(m,res_norm,n);
%     if sum(s)<= res_norm
%     l=l+1;
%     res_vec(l,1:m)=s;
%     end
%     n=n+1;
% end
% 
% res_size=size(res_vec);
% 
% while n32<= T
%      g=randi([1,res_size(1)]);
%      D(1:m,n32)=res_vec(g,1:m);
%      n32=n32+1;
%  end
%load New_D_N6_m3.mat
% 
% while n32<= T
%     g=1*rand(m,1);
%     g=round(g);
%     if sum(g)<=N
%         D(1:m,n32)=g;
%         n32=n32+1;
%     end
% end

%D=[0 1 0;1 0 1];

%%

%C=[0.9 1 1.1];
%n1=1;
%y2=Addition(D(1:m,T)',C);
% while n1<= NS(1,1)
%     [ok] = state_difference(S(n1,1:end),D(1:m,end));
%     if ok==1 
%         y1=Addition(S(n1,1:end),C);
%         Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end), y1];
%     else
%         Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end), 0];
%     end
%     n1=n1+1;
% end


%D=[1 0 0 1 0 1;0 1 1 1 0 0;0 1 0 1 1 1];
% C=[0.9 1 1.1];
% n1=1;
% while n1<= NS(1,1)
%     [ok] = state_difference(S(n1,1:end),D(:,end));
%     if ok==1
%        Addition(S(n1,1:end),C);
%        Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end) Addition(S(n1,1:end),C)];
%     else
%         Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end) 0];
%     end
%     %Addition(S(n1,1:end),C);
%     %Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end) Addition(S(n1,1:end),C)];
%     n1=n1+1;
% end


%% 
a=1;
n1=1;
n2=1;
k=1;
%k=T-1;

loop_indicator=0;

%%
%while k>= 1
    %lambda=[0.5+0.1*cos(k*pi/3);0.3+.05*sin(k*pi/7)];
    %mu=[0.3-0.1*cos(k*pi/3);0.4-0.05*sin(k*pi/7)];
    k;
    %D_k=D(:,k);
    %D_k_p=D(:,k+1);
    %%
    if norm(D_k_p+F_k_p,1)==N & norm(D_k+F_k,1)<N
   
        if sum(state) <= N
           transition=D_k_p'+F_k_p';
           state=transition;
        end
       
    end
    
    %%
    if norm(D_k_p+F_k_p,1)==N & norm(D_k+F_k,1)==N
        
         if state == D_k'+F_k'
                transition=D_k_p'+F_k_p';
                state=transition;
         end
           
    end
    
    
    %%
    
    
    if norm(D_k+F_k,1)==N & norm(D_k_p+F_k_p,1)< N & equal_less==0
        
        equal_less=1;
        
        n1;
        state;
        
        if Decision~=0 & Decision < m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
        end
        
        if Decision == m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
        end
        
        
        if Decision == m+2
                ok_2=0;
                er=0;
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=state-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state_unbounded(SS_p,D_k_p,F_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        state = shifted_transition_unbounded(SS_p,D_k_p,F_k_p,m);
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end 
                    
                    if er==m & ok_2==0
                        ok_2=1;
                        transition_6 = shifted_transition_unbounded(zeros(1,m),D_k_p,F_k_p,m);
                        state = transition_6;
                    end
                    
                end
        end
            
    end
    %%
    
    
if norm(D_k_p+F_k_p,1) < N & norm(D_k+F_k,1) < N & less_less==0 & equal_less==0
        
        less_less=1;
        state;
            
            if Decision~=0 & Decision < m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
            if Decision == m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
            if Decision == m+2
                ok_2=0;
                er=0;
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=state-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state_unbounded(SS_p,D_k_p,F_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        state = shifted_transition_unbounded(SS_p,D_k_p,F_k_p,m);
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        ok_2=1;
                        transition_6 = shifted_transition_unbounded(zeros(1,m),D_k_p,F_k_p,m);
                        state = transition_6;
                    end
                    
                end
            end
            
        end
         
        
        %%
        
        if norm(D_k+F_k,1)==0 & less_less==0 & equal_less==0 & curr_demand==0
            
            curr_demand=1;
            
            state;
            D_k';
            D_k_p';
            
            
            if Decision~=0 & Decision < m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
            if Decision == m+1
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
            if Decision == m+2
                ok_2=0;
                er=0;
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=state-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state_unbounded(SS_p,D_k_p,F_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        state = shifted_transition_unbounded(SS_p,D_k_p,F_k_p,m);
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                   if er==m & ok_2==0
                        ok_2=1;
                        transition_6 = shifted_transition_unbounded(zeros(1,m),D_k_p,F_k_p,m);
                        state = transition_6;
                    end 
                    
                end
            end
            
        end

        %%
        
        
        
        if (sum(state)==N) & less_less==0 & equal_less==0 & curr_demand==0
            
            full_system=1;
            
            D_k;
            D_k_p;
            
            state;
            
            if norm(D_k+F_k,1) == 0
                
                
                
                if Decision == m+1
                S_p=S(index_state,1:end)-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);

                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
                end
                
                 if Decision == m+2
                ok_2=0;
                er=0;
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=state-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state_unbounded(SS_p,D_k_p,F_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        state = shifted_transition_unbounded(SS_p,D_k_p,F_k_p,m);
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end                    
                end
            end
            end
            
            if norm(D_k+F_k,1) ~= 0
               state;
               
                if Decision~=0 & Decision < (m+1)
                    S_p=state-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                    state_path=path_creation(transition,combined_probability,f);
                    state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
                end
                if Decision == (m+1)
                    S_p=state-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                    state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                    state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
                end
                
                if Decision == m+2
                ok_2=0;
                er=0;
                S_p=state-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=state-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state_unbounded(SS_p,D_k_p,F_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        state = shifted_transition_unbounded(SS_p,D_k_p,F_k_p,m);
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        ok_2=1;
                        transition_6 = shifted_transition_unbounded(zeros(1,m),D_k_p,F_k_p,m);
                        state = transition_6;
                    end 
                    
                end
            end
                
            end
        end
        
        next_state=state;
        
        
        
        
 
        
   
    
   
    
%end




