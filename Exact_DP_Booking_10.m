function [Revenue S] = Exact_DP_Booking_10(N,m,T,D,lambda_max,lambda_min,mu_max,mu_min,C)


%%================ Exact DP Booking Easter Period 

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

n=0;
l=0;
%m=3;
%N=4;
while n<=((N+1)^m - 1)
    s=base(m,N,n);
    if sum(s)<= N
    l=l+1;
    S(l,1:m)=s;
    end
    n=n+1;
end

% NS: size of entire state space
NS=size(S);


%% Generating the reservation vector
% 
n32=1;
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
n1=1;
%y2=Addition(D(1:m,T)',C);
while n1<= NS(1,1)
    [ok] = state_difference(S(n1,1:end),D(1:m,end));
    if ok==1 
        y1=Addition(S(n1,1:end),C);
        Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end), y1];
    else
        Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end), 0];
    end
    n1=n1+1;
end


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


%% In this subsection we apply exact DP for each state and the entire horizon
% we use variable "a" as a temporarly variable for price
% "k" is the time slot counter
% "Revenue" is the variable that saves the value function for each time
% slot, it is a 3 dimension matrix; also it saves the optimal decision 
% for each state

a=1;
n1=1;
n2=1;

k=T-1;
Revenue(1:NS(1,1),1:(m+1),T)=Terminal_revenue(1:NS(1,1),1:(m+1));
Revenue(1:NS(1,1),(m+2),T)=0;

% This subsection is the main loop
% we check whether or not each price is zero 
% loop on the variable "k" represents the horizon which is the backward
% procedure; loop on variable "n1" represents the state space that evaluate
% the policy for each state; loop on variable "a" represents evaluating a
% price for each state
%=====================
% stateanalysis(s,N,m,a): It is a function whose inputs is the state "s",
% the capacity "N", number of prices "m", and single price "a=c_i"; whose
% output is possible transitions from a state if an action has taken; the
% other output is the probability of going from a state to other state by 
% taking an action, 
% e.g. [transition combined_probability]=stateanalysis([1 1],3,2,1)
%=====================
% at each time instant "k" and each state, with having possible transition
% and the probability of each transition, we want to know the value
% function at the previous step! so we have written a function named
% "pair_finding" whose inputs are "transition", "Revenue" of the previous
% step, "size of state space NS(1,1)", "probability of each transition"

% possible_Revenue(1,a): It is a vector that saves the value fucntion for
% taking a price "a=c_i" its size is 1*m i.e. at each state we evaluate "m"
% different policy and save it into this variable, the DP chooses the best
% policy from this vector for each state at each time instant and inserts
% into the "Revenue"

loop_indicator=0;

while k>= 1
    %lambda=[0.5+0.1*cos(k*pi/3);0.3+.05*sin(k*pi/7)];
    %mu=[0.3-0.1*cos(k*pi/3);0.4-0.05*sin(k*pi/7)];
    k;
    D_k=D(:,k);
    D_k_p=D(:,k+1);
    %%
    if norm(D_k_p,1)==N & norm(D_k,1)<N
        while n1<= NS(1,1)
            [ok] = state_difference(S(n1,1:end),D_k);
            if ok == 1
                if sum(S(n1,1:end))<=N
                    transition=D_k_p';
                    combined_probability=ones(1,m);
                    x=Revenue(1:end,1:end,k+1);
                    Rev=Pair_finding_2(transition,x,NS(1,1),combined_probability);
                    Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                    Revenue(n1,1:m,k)=S(n1,1:end);
                    Revenue(n1,(m+2),k)=0;
                end
            else
                    Revenue(n1,1:m,k)=S(n1,1:end);
                    Revenue(n1,(m+1),k)=0;
                    Revenue(n1,(m+2),k)=0;
            end
            n1=n1+1;
        end
        k=k-1;
        n1=1;
    end
    
    %%
    if norm(D_k_p,1)==N & norm(D_k,1)==N
        while n1<= NS(1,1)
            
            if (S(n1,1:end)) == D_k'
                transition=D_k_p';
                combined_probability=ones(1,m);
                x=Revenue(1:end,1:end,k+1);
                Rev=Pair_finding_2(transition,x,NS(1,1),combined_probability);
                Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=0;
            else
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+1),k)=0;
                Revenue(n1,(m+2),k)=0;
            end
            n1=n1+1;
        end
        k=k-1;
        n1=1;
    end
    
    
    %%
    
    if norm(D_k,1)==N & norm(D_k_p,1)< N
        while n1<= NS(1,1)
            if S(n1,1:end)== D_k'
                while a <= m
                    if C(a)~=0
                        S(n1,1:end);
                        [ok] = state_difference(S(n1,1:end),D_k);
                        
                       if ok==1
                        S_p=S(n1,1:end)-D_k';
                        [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                        [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                        [state_status]=eligible_state(transition,D_k_p,m,N);
                        x=Revenue(1:end,1:end,k+1);
                        Revenue(1:end,1:end,k+1);
                       
                       if state_status==1
                           transition_2 = shifted_transition(transition,D_k_p,m);
                           Rev=Pair_finding(transition_2,x,NS(1,1),combined_probability);
                           possible_Revenue(1,a)=Rev;
                        else
                           possible_Revenue(1,a)=0;
                       end
                       
                        else
                        possible_Revenue(1,a)=0;
                    end
                        
                    end
                    a=a+1;
                    indicator_rej=1;
                end
                
            if indicator_rej==1
                [ok] = state_difference(S(n1,1:end),D_k);
                if ok==1
                    S_p=S(n1,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    x_re=Revenue(1:end,1:end,k+1);
                    if state_status==1
                           transition_3 = shifted_transition(transition_re,D_k_p,m);
                           Rev=Pair_finding(transition_3,x,NS(1,1),combined_probability);
                           possible_Revenue(1,1+m)=Rev;
                       else
                           possible_Revenue(1,1+m)=0;
                    end
                else
                    possible_Revenue(1,1+m)=0;
                end
                
            end
            [ok] = state_difference(S(n1,1:end),D_k);
            %norm(S(n1,1:end),1)>= D_k
            if ok==1
                indicator_rej=0;
                [max_revenue_state price]=max(possible_Revenue);
                Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+ max_revenue_state;
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=price;
            else
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=0;
            end
            
            end
            n1=n1+1;
            a=1;
            indicator_rej=0;
        end
        k=k-1;
        n1=1;
    end
    
    
    %%
    if norm(D_k_p,1) < N & norm(D_k,1) < N
    while n1<= NS(1,1) 
        
        if norm(D_k,1)> 0
            k;
        if sum(S(n1,1:end))<N
            state_feasibility_3=0;
            while a <= (m) 
                if C(a)~=0
                    S(n1,1:end);
                    %D_k=D(:,k);
                    %D_k_p=D(:,k+1);
                    [ok] = state_difference(S(n1,1:end),D_k);
                    
                    if ok==1
                       S_p=S(n1,1:end)-D_k';
                       [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);

                       [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition,D_k_p,m,N);
                       x=Revenue(1:end,1:end,k+1);
                       Revenue(1:end,1:end,k+1);
                       
                       if state_status==1
                           state_feasibility_3=state_feasibility_3+1;
                           transition_2 = shifted_transition(transition,D_k_p,m);
                           Rev=Pair_finding(transition_2,x,NS(1,1),combined_probability);
                           possible_Revenue(1,a)=Rev;
                       else
                           possible_Revenue(1,a)=0;
                       end
                       
                       %Rev=Pair_finding(transition,x,NS(1,1),combined_probability);
                       %possible_Revenue(1,a)=Rev;
                       
                    else
                        possible_Revenue(1,a)=0;
                    end
                   
                end
                a=a+1;
                indicator_rej=1;
            end
            
            if indicator_rej==1
                [ok] = state_difference(S(n1,1:end),D_k);
                if ok==1
                    S_p=S(n1,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    x_re=Revenue(1:end,1:end,k+1);
                    if state_status==1
                           state_feasibility_3=state_feasibility_3+1;
                           transition_3 = shifted_transition(transition_re,D_k_p,m);
                           Rev=Pair_finding(transition_3,x,NS(1,1),combined_probability);
                           possible_Revenue(1,1+m)=Rev;
                       else
                           possible_Revenue(1,1+m)=0;
                    end
                else
                    possible_Revenue(1,1+m)=0;
                end
                
            end
            [ok] = state_difference(S(n1,1:end),D_k);
            %norm(S(n1,1:end),1)>= D_k
            if ok==1 & state_feasibility_3~=0
                indicator_rej=0;
                [max_revenue_state price]=max(possible_Revenue);
                Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+ max_revenue_state;
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=price;
            end
            
            if ok==1 & state_feasibility_3==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=S(n1,1:end)-D_k';
                %[min_variation_rate min_price_rate]=min(variation_rate);
                %[min_variation_rate
                %min_price_rate]=min(variation_rate_matrix)
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    x=Revenue(1:end,1:end,k+1);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    end
                    
                end
            end
            
            %else
                %Revenue(n1,1:m,k)=S(n1,1:end);
                %Revenue(n1,(m+2),k)=0;
            %end
            
        end
        
        end
        
        %%
        
        if norm(D_k,1)==0
            state_feasibility=0;
            k;
        if sum(S(n1,1:end))<N
            while a <= (m) 
                if C(a)~=0
                    S(n1,1:end);
                    %D_k=D(:,k);
                    %D_k_p=D(:,k+1);
                    k;
                    [ok] = state_difference(S(n1,1:end),D_k);
                    
                    if ok==1
                       S_p=S(n1,1:end)-D_k';
                       %S_p=[2;0;0];
                       %a=1
                       [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                       [transition_zero combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition_zero,D_k_p,m,N);
                       
                       k;
                       x=Revenue(1:end,1:end,k+1);
                       Revenue(1:end,1:end,k+1);
                       
                       if state_status==1
                           S_p;
                           k;
                           state_feasibility=state_feasibility+1;
                           transition_2_zero = shifted_transition(transition_zero,D_k_p,m);
                           Rev=Pair_finding(transition_2_zero,x,NS(1,1),combined_probability);
                           possible_Revenue(1,a)=Rev;
                       else
                           possible_Revenue(1,a)=0;
                       end
                       
                    else
                        possible_Revenue(1,a)=0;
                    end
                   
                end
                a=a+1;
                indicator_rej=1;
            end
            
            if indicator_rej==1
                [ok] = state_difference(S(n1,1:end),D_k);
                if ok==1
                    S_p=S(n1,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_re_zero combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re_zero,D_k_p,m,N);
                    x_re=Revenue(1:end,1:end,k+1);
                    if state_status==1
                           state_feasibility=state_feasibility+1;
                           S_p;
                           k;
                           transition_3_zero = shifted_transition(transition_re_zero,D_k_p,m);
                           Rev=Pair_finding(transition_3_zero,x,NS(1,1),combined_probability);
                           possible_Revenue(1,1+m)=Rev;
                       else
                           possible_Revenue(1,1+m)=0;
                    end
                else
                    possible_Revenue(1,1+m)=0;
                end
                
            end
            
            % here we have to define new variable to be aware of revenue
            % function
            %state_feasibility
            
            [ok] = state_difference(S(n1,1:end),D_k);
            %norm(S(n1,1:end),1)>= D_k
            if ok==1 & state_feasibility~=0
                S_p;
                possible_Revenue;
                indicator_rej=0;
                [max_revenue_state price]=max(possible_Revenue);
                k;
                
                Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+ max_revenue_state;
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=price;
            end
            %else
             %   indicator_rej=0;
             %  Revenue(n1,1:m,k)=S(n1,1:end);
             %   Revenue(n1,(m+2),k)=0;
             
             if ok==1 & state_feasibility==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=S(n1,1:end)-D_k';
                %[min_variation_rate min_price_rate]=min(variation_rate);
                %[min_variation_rate
                %min_price_rate]=min(variation_rate_matrix)
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    x=Revenue(1:end,1:end,k+1);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    end
                    
                    
                end
            end
                         
            end
            
        end
        
           %state_feasibility=0; 
        %end
        
        
        
        %%
        
        if sum(S(n1,1:end))==N
            %[ok] = state_difference(S(n1,1:end),D_k);
            %norm(S(n1,1:end),1)>= D_k
            if norm(D_k,1) == 0
                state_feasibility_2=0;
                [ok] = state_difference(S(n1,1:end),D_k);
            if ok==1
                S_p=S(n1,1:end)-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                [state_status]=eligible_state(transition_ca,D_k_p,m,N);
                x_ca=Revenue(1:end,1:end,k+1);
                transition_ca;
                if state_status==1
                           transition_4 = shifted_transition(transition_ca,D_k_p,m);
                           state_feasibility_2=state_feasibility_2+1;
                           Rev_ca=Pair_finding(transition_4,x_ca,NS(1,1),combined_probability_ca);
                           possible_Revenue_ca=Rev_ca;
                           Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+ possible_Revenue_ca;
                           Revenue(n1,1:m,k)=S(n1,1:end);
                           Revenue(n1,(m+2),k)=m+1;
                %else
                %        possible_Revenue(1,1+m)=0;
                %         Revenue(n1,1:m,k)=S(n1,1:end);
                %         Revenue(n1,(m+2),k)=0;
                         
                end
                
                %Rev_ca=Pair_finding(transition_ca,x_ca,NS(1,1),combined_probability_ca);
                %possible_Revenue_ca=Rev_ca;
                %Revenue(n1,(m+1),k)=Terminal_revenue(n1,(m+1))+ possible_Revenue_ca;
                %Revenue(n1,1:m,k)=S(n1,1:end);
                %Revenue(n1,(m+2),k)=m+1;
            %else
            %    Revenue(n1,1:m,k)=S(n1,1:end);
            %    Revenue(n1,(m+2),k)=0;
            
            if state_feasibility_2==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=S(n1,1:end)-D_k';
                %[min_variation_rate min_price_rate]=min(variation_rate);
                %[min_variation_rate
                %min_price_rate]=min(variation_rate_matrix)
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    x=Revenue(1:end,1:end,k+1);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    end
                    
                end
            end
            
            
            
            
            
            end
            
            end
            
            if norm(D_k,1) ~= 0
                %[ok] = state_difference(S(n1,1:end),D_k);
                state_feasibility_2=0;
                a=1;
                while a <= (m) 
                if C(a)~=0
                    S(n1,1:end);
                    %D_k=D(:,k);
                    %D_k_p=D(:,k+1);
                    [ok] = state_difference(S(n1,1:end),D_k);
                    
                    if ok==1
                       S_p=S(n1,1:end)-D_k';
                       [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                       [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition,D_k_p,m,N);
                       x=Revenue(1:end,1:end,k+1);
                       Revenue(1:end,1:end,k+1);
                       
                       if state_status==1
                           state_feasibility_2=state_feasibility_2+1;
                           transition_2 = shifted_transition(transition,D_k_p,m);
                           Rev=Pair_finding(transition_2,x,NS(1,1),combined_probability);
                           possible_Revenue(1,a)=Rev;
                       else
                           possible_Revenue(1,a)=0;
                       end
                       
                    else
                        possible_Revenue(1,a)=0;
                    end
                   
                end
                a=a+1;
                indicator_rej=1;
            end
                
                if indicator_rej==1
                [ok] = state_difference(S(n1,1:end),D_k);
                if ok==1
                    S_p=S(n1,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    
                    x_re=Revenue(1:end,1:end,k+1);
                    if state_status==1
                            state_feasibility_2=state_feasibility_2+1;
                           transition_3 = shifted_transition(transition_re,D_k_p,m);
                           Rev=Pair_finding(transition_3,x,NS(1,1),combined_probability);
                           possible_Revenue(1,1+m)=Rev;
                       else
                           possible_Revenue(1,1+m)=0;
                    end
                 else
                    possible_Revenue(1,1+m)=0;
                end
                end
                
                [ok] = state_difference(S(n1,1:end),D_k);
            %norm(S(n1,1:end),1)>= D_k
            if ok==1 & state_feasibility_2~=0
                indicator_rej=0;
                [max_revenue_state price]=max(possible_Revenue);
                Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+ max_revenue_state;
                Revenue(n1,1:m,k)=S(n1,1:end);
                Revenue(n1,(m+2),k)=price;
            %else
            %    Revenue(n1,1:m,k)=S(n1,1:end);
            %    Revenue(n1,(m+2),k)=0;
            end
              
            if ok==1 & state_feasibility_2==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=S(n1,1:end)-D_k';
                %[min_variation_rate min_price_rate]=min(variation_rate);
                %[min_variation_rate
                %min_price_rate]=min(variation_rate_matrix)
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    x=Revenue(1:end,1:end,k+1);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    else
                        er=er+1;
                        variation_rate_matrix(min_price_rate,:)=[];
                    end
                    if er==m & ok_2==0
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        Rev=Pair_finding(transition_6,x,NS(1,1),combined_probability);
                        Revenue(n1,(m+1),k)=Addition(S(n1,1:end),C)+Rev;
                        Revenue(n1,1:m,k)=S(n1,1:end);
                        Revenue(n1,(m+2),k)=m+2;
                    end
                end
            end
            
            
            
            end
            
        end
        n1=n1+1;
        a=1;
        indicator_rej=0;
    end
    k=k-1;
    n1=1;
    end
    
end




