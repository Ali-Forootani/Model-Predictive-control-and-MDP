
%%================ Exact DP Booking Easter Period 

clear all
clc
tic
%%
m=input('Enter the number of prices:')
N=input('Enter the parking capacity:')



%% Time Horizon, price vector,
% T: variable for the horizon
%  C: Price vector

%D_k=[0;0;0];
%D_k_p=[1;1;1];
%%
%lambda=[0.6;0.5;0.3];
%mu=[0.2;0.2;0.4];
lambda=[0.6;0.5;0.3;0.2];
mu=[0.2;0.2;0.4;0.4];
%%
%lambda=[0.6;0.5];
%mu=[0.2;0.2];
T=30;
%C=[0.9 1];
%C=[0.9 1 1.1];
C=[0.9 1 1.1 1.2];
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
%==========================================
% while n<=((N+1)^m - 1)
%     s=base(m,N,n);
%     if sum(s)<= N
%     l=l+1;
%     S(l,1:m)=s;
%     end
%     n=n+1;
% end

% NS: size of entire state space
% NS=size(S);


%% Generating the reservation vector
% 
n32=1;
% 
% while n32<= T
%     g=1*rand(m,1);
%     g=round(g);
%     if sum(g)<=N
%         D(1:m,n32)=g;
%         n32=n32+1;
%     end
% end

%load D_N6_m3.mat

%D=[0 1 0;1 0 1];

%%

%C=[0.9 1 1.1];
n1=1;


%% In this subsection we apply exact DP for each state and the entire horizon
% we use variable "a" as a temporarly variable for price
% "k" is the time slot counter
% "Revenue" is the variable that saves the value function for each time
% slot, it is a 3 dimension matrix; also it saves the optimal decision 
% for each state

a=1;
n1=1;
n2=1;

%k=T-1;

%Revenue(1:NS(1,1),1:(m+1),T)=Terminal_revenue(1:NS(1,1),1:(m+1));
%Revenue(1:NS(1,1),(m+2),T)=0;

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
initial_state=[8 9 10 10];
dis_factor=0.996;
horizon=30;
load('E:\Education\park_pricing\1.16.2017\Reservation_examples\ADP_reservation_example\r_hat_N50_m4.mat')
load('E:\Education\park_pricing\1.16.2017\Reservation_examples\ADP_reservation_example\D_N50_m4.mat')

k=7;
l=21;

Decision=[];
price_vector=[0.9 1 1.1 1.2];
s=initial_state;
kk=1;
l_m=3;
n=1;
Expected_Revenue(n)=0;

num_exper=10000;

%%
while n <num_exper
    initial_state=[1 3 1 9]+D(:,k)';
    s=initial_state;
    n1=k;
    Expected_Revenue(n)=0;
while n1 <= l-2
    %lambda=[0.5+0.1*cos(k*pi/3);0.3+.05*sin(k*pi/7)];
    %mu=[0.3-0.1*cos(k*pi/3);0.4-0.05*sin(k*pi/7)];
    %%
    k;
    D_k=D(:,n1);
    D_k_p=D(:,n1+1);
    f=rand;
    Expected_Revenue(n)=Expected_Revenue(n)+Addition(s,C);
    s;
    kk=kk+1;
    
    %%
    if norm(D_k_p,1)==N & norm(D_k,1)<N
            [ok] = state_difference(s,D_k);
            if ok == 1
                if sum(s)<=N
                    %%
                    transition=D_k_p';
                    next_state=D_k_p';
                    Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                end
            else
            end
        n1=n1+1;
        s=next_state;
    end
    
    %%
    if norm(D_k_p,1)==N & norm(D_k,1)==N
            if (s) == D_k'                
                %%
                    transition=D_k_p';
                    next_state=D_k_p';
                    Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
            else
            end
        n1=n1+1;
        s=next_state;
    end    
    %%
    if norm(D_k,1)==N & norm(D_k_p,1)< N
            if s== D_k'
                while a <= m
                    if C(a)~=0
                        [ok] = state_difference(s,D_k);
                       if ok==1
                        S_p=s-D_k';
                        [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                        [state_status]=eligible_state(transition,D_k_p,m,N);
                       if state_status==1
                           transition_2 = shifted_transition(transition,D_k_p,m);
                           %%
                           [next_state]=path_creation(transition_2,combined_probability,f,a);
                           possible_nextstate(a,:)=next_state;
                           
                           
                           %%
                           
                           ph_kprime=phi_cal2(m,next_state);
                           
                           %%
                           possible_revenue(a)=prod(ph_kprime,r_hat);
                           
                           %%
                           g_state=g_cal(s,price_vector);
                            
                       else
                           possible_revenue(a)=0;
                           
                       end
                       
                       else
                           possible_revenue(a)=0;
                           
                    end
                        
                    end
                    a=a+1;
                    indicator_rej=1;
                end
                
            if indicator_rej==1
                [ok] = state_difference(S(n1,1:end),D_k);
                if ok==1
                    S_p=s-D_k';
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    if state_status==1
                           transition_3 = shifted_transition(transition_re,D_k_p,m);
                           
                           %%
                           [next_state]=path_creation_re(transition_3,combined_probability_re,f);
                           possible_nextstate(m+1,:)=next_state;
                           
                           %%
                           ph_kprime=phi_cal2(m,next_state);
                           
                           %%
                           possible_revenue(m+1)=prod(ph_kprime,r_hat);
                           
                    else
                           possible_revenue(m+1)=0;
                           
                    end
                else
                    possible_revenue(m+1)=0;
                end
            end
            [ok] = state_difference(s,D_k);
            if ok==1
                indicator_rej=0;
                %r_norm_path(k,1:m,l_m)=r_norm;
                [r_max_instant Decision]=max(possible_revenue);
                Decision;
                s;
                next_state=possible_nextstate(Decision,:);
                Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
            else
            end
            end
            n1=n1+1;
            a=1;
            indicator_rej=0;
        s=next_state;
    end
    
    %%
    if norm(D_k_p,1) < N & norm(D_k,1) < N
        
        if norm(D_k,1)> 0
            k;
        if sum(s)<N
            state_feasibility_3=0;
            while a <= (m) 
                if C(a)~=0
                    [ok] = state_difference(s,D_k);
                    if ok==1
                       S_p=s-D_k';
                       [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition,D_k_p,m,N);
                       if state_status==1
                           state_feasibility_3=state_feasibility_3+1;
                           transition_2 = shifted_transition(transition,D_k_p,m);
                           %%
                           [next_state]=path_creation(transition_2,combined_probability,f,a);
                           possible_nextstate(a,:)=next_state;
                           %%
                           ph_kprime=phi_cal2(m,next_state);
                           %%
                           possible_revenue(a)=dot(ph_kprime,r_hat);
                       else
                           possible_revenue(a)=0;
                       end
                    else
                        possible_revenue(a)=0;
                    end
                end
                a=a+1;
                indicator_rej=1;
            end
            
            if indicator_rej==1
                [ok] = state_difference(s,D_k);
                if ok==1
                    S_p=s-D_k';
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    if state_status==1
                           
                           %%
                           state_feasibility_3=state_feasibility_3+1;
                           transition_3 = shifted_transition(transition_re,D_k_p,m);
                           
                           %%
                           [next_state]=path_creation_re(transition_3,combined_probability_re,f);
                           possible_nextstate(m+1,:)=next_state;
                           
                           %%
                           ph_kprime=phi_cal2(m,next_state);
                           possible_revenue(m+1)=dot(ph_kprime,r_hat);
                           
                    else
                           possible_revenue(m+1)=0;
                          
                    end
                else
                           possible_revenue(m+1)=0;
                           
                end
                
            end
            [ok] = state_difference(s,D_k);
            if ok==1 & state_feasibility_3~=0
                indicator_rej=0;
                [r_max_instant Decision]=max(possible_revenue);
                Decision;
                s;
                next_state=possible_nextstate(Decision,:);
                s=next_state;
                Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
               
            end
            
            if ok==1 & state_feasibility_3==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=s-D_k';
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                       
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(m+2)=dot(ph_kprime,r_hat);
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);    
                        s=next_state;
                        
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(m+2)=dot(ph_kprime,r_hat);
                        
                        %%
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        s=next_state;
                       
                    end
                    
                end
            end
            
        end
        
        end
        
        %%
        
        if norm(D_k,1)==0
            state_feasibility=0;
            k;  
        if sum(s)<N
            while a <= (m) 
                if C(a)~=0
                    s;
                    k;
                    [ok] = state_difference(s,D_k);
                    if ok==1
                       S_p=s-D_k';
                       [transition_zero combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition_zero,D_k_p,m,N);
                       k;
                       
                       if state_status==1
                        S_p;
                        k;
                        state_feasibility=state_feasibility+1;
                        transition_2_zero = shifted_transition(transition_zero,D_k_p,m);
                           
                        %%
                        [next_state]=path_creation(transition_zero,combined_probability,f,a);
                        possible_nextstate(a,:)=next_state;
                           
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(a)=dot(ph_kprime,r_hat);
                        
                       else
                          possible_Revenue(a)=0;
                       end
                       
                    else
                        possible_Revenue(a)=0;
                    end
                   
                end
                a=a+1;
                indicator_rej=1;
            end
            
            if indicator_rej==1
                [ok] = state_difference(s,D_k);
                if ok==1
                    S_p=s-D_k';
                    [transition_re_zero combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re_zero,D_k_p,m,N);
                    if state_status==1
                        state_feasibility=state_feasibility+1;
                        S_p;
                        k;
                        transition_3_zero = shifted_transition(transition_re_zero,D_k_p,m);
                        %%
                        [next_state_ca]=path_creation_re(transition_3_zero,combined_probability_re,f);   
                        possible_nextstate(m+1,:)=next_state_ca;
                           
                        %%
                        ph_kprime=phi_cal2(m,next_state_ca);
                        
                        %%
                        
                        possible_revenue(m+1)=dot(ph_kprime,r_hat);
                        
                    else
                        possible_revenue(m+1)=0;
                    end
                else
                    possible_revenue(m+1)=0;
                end
                
            end
            [ok] = state_difference(s,D_k);
            if ok==1 & state_feasibility~=0
            S_p;
            current_state=s;
            [r_max_instant Decision]=max(possible_revenue);   
            next_state=possible_nextstate(Decision,:);   
            Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
            s=next_state;    
            end
             if ok==1 & state_feasibility==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=s-D_k';
                SS_p=S_p;
                while ok_2==0
                    %SS_p=S_p;
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(m+2)=prod(ph_kprime,r_hat);
                        s=next_state;    
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        %%
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(m+2)=dot(ph_kprime,r_hat);
                        s=next_state;
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        
                    end
                    
                end
            end
             
            end
            
        end
        
        %%
        if sum(s)==N
            
            if norm(D_k,1) == 0
                state_feasibility_2=0;
                [ok] = state_difference(s,D_k);
            if ok==1
                S_p=s-D_k';
                [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                [state_status]=eligible_state(transition_ca,D_k_p,m,N);
                if state_status==1
                        %%
                        transition_4 = shifted_transition(transition_ca,D_k_p,m);
                        state_feasibility_2=state_feasibility_2+1;
                        [next_state_ca]=path_creation_re(transition_4,combined_probability_ca,f);   
                        possible_nextstate(m+1,:)=next_state_ca;
                           
                        %%
                        ph_kprime=phi_cal2(m,next_state_ca);
                        
                        %%
                        possible_revenue(m+1)=dot(ph_kprime,r_hat);
                        
                end
            
            if state_feasibility_2==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=s-D_k';
                
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(m+2)=dot(ph_kprime,r_hat);
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        s=next_state;                            
                    else
                        variation_rate_matrix(min_price_rate,:)=[];
                        er=er+1;
                    end
                    
                    if er==m & ok_2==0
                        %%
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        
                        %%
                        [next_state]=transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        s=next_state;    
                      
                    end
                    
                end
            end
        end
     end
            
            if norm(D_k,1) ~= 0
                state_feasibility_2=0;
                a=1;
                while a <= (m) 
                if C(a)~=0
                    
                    [ok] = state_difference(s,D_k);
                    if ok==1
                       S_p=s-D_k';
                       [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
                       [state_status]=eligible_state(transition,D_k_p,m,N);
                       
                       if state_status==1
                        %%
                        state_feasibility_2=state_feasibility_2+1;
                        transition_2 = shifted_transition(transition,D_k_p,m);
                         
                        %%
                        [next_state]=path_creation(transition_2,combined_probability,f,a);
                        possible_nextstate(a,:)=next_state;
                           
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        
                        %%
                        possible_revenue(a)=dot(ph_kprime,r_hat);
                        
                       
                       else
                           possible_revenue(a)=0;
                          
                       end
                       
                    else
                        possible_revenue(a)=0;
                       
                    end
                   
                end
                a=a+1;
                indicator_rej=1;
            end
                
                if indicator_rej==1
                [ok] = state_difference(s,D_k);
                if ok==1
                    S_p=s-D_k';
                    [transition_re combined_probability_re]=stateanalysis_re_bo_ver_3(S_p,N,m,mu,D_k,D_k_p);
                    [state_status]=eligible_state(transition_re,D_k_p,m,N);
                    
                    if state_status==1
                        %%   
                        state_feasibility_2=state_feasibility_2+1;
                        transition_3 = shifted_transition(transition_re,D_k_p,m);
                        
                        %%
                        [next_state_re]=path_creation_re(transition_3,combined_probability_re,f);   
                        possible_nextstate(m+1,:)=next_state_re;
                           
                           %%
                           
                           ph_kprime=phi_cal2(m,next_state_re);
                           
                           %%
                           possible_revenue(m+1)=dot(ph_kprime,r_hat);
                           
                    else
                          possible_revenue(m+1)=0;
                          
                    end
                else
                    
                    possible_revenue(m+1)=0;
                    
                end
            end
                
            [ok] = state_difference(s,D_k);
            
            if ok==1 & state_feasibility_2~=0
                indicator_rej=0;
%                 for i=1:m+1
%                     instant_rev(i)=g_state+dis_factor*phi_cal2(m,possible_nextstate(i,:))'*r(:,i);
%                 end
       
            current_state=s;
            [r_max_instant Decision]=max(possible_revenue);
            next_state=possible_nextstate(Decision,:);   
            s=next_state;
            Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
            
            end
              
            if ok==1 & state_feasibility_2==0
                indicator_rej=0;
                ok_2=0;
                er=0;
                variation_rate=lambda-mu;
                for n35=1:m
                    variation_rate_matrix(n35,:)=[variation_rate(n35),n35];
                end
                S_p=s-D_k';
                SS_p=S_p;
                while ok_2==0
                    [min_variation_rate min_price_rate]=min(variation_rate_matrix(:,1));
                    SS_p(1,min_price_rate)=0;
                    [state_status_5]=eligible_state(SS_p,D_k_p,m,N);
                    if state_status_5==1
                        ok_2=1;
                        %%
                        transition_6 = shifted_transition(SS_p,D_k_p,m);
                        combined_probability_6=ones(1,m);
                        next_state = transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        s=next_state;
                    
                    else
                        er=er+1;
                        variation_rate_matrix(min_price_rate,:)=[];
                    end
                    if er==m & ok_2==0
                        %%
                        transition_6 = shifted_transition(zeros(1,m),D_k_p,m);
                        combined_probability_6=ones(1,m);
                        next_state = transition_6;
                        
                        %%
                        ph_kprime=phi_cal2(m,next_state);
                           
                        %%
                        Expected_Revenue(n)=Expected_Revenue(n)+Addition(next_state,C);
                        s=next_state;
                        
                    end
                end
            end
            
            end
            
        end
        %n1=n1+1;
        %a=1;
        %indicator_rej=0;
    %end
    %k=k-1;
    %n1=1;
    %n1=n1+1;
    %decision=[];
    %possible_nextstate=[];
    end
    n1=n1+1;
    a=1;
    indicator_rej=0;
    decision=[];
    possible_nextstate=[];
    possible_revenue=[];
end
n1=k;
%kk=0;
state_visit(n,1:m)=s;
n=n+1;
end

%% mean expected revenue for the time duration of k-(l-2)
%num_exper=499;
mean_expected_revenue=sum(Expected_Revenue)/num_exper;

%% Computing the probability distribution for the time slot "l-1"

st_size=size(state_visit);
n2=1;
%num_exper=499;
n4=1;
n3=2;
state_space(n4,1:m)=state_visit(1,1:m);
while n3 <= st_size(1)-2
    s=state_visit(n3,1:m);
    rep_ok=Pair_finding_ADP(state_space,s);
    if rep_ok==0
        n4=n4+1;
        state_space(n4,1:m)=s;
    end
    n3=n3+1;
end

%% Probability Distribution at the time slot "l-1"
state_space_size=size(state_space);
n2=1;
while n2 <= state_space_size(1)
    D_k = D(:,l-1);
    [ok] = state_difference(state_space(n2,1:m),D_k);
    if ok==1
        [frequency_visit prob_visit]=Pair_finding_dis(state_visit,state_space(n2,1:m),num_exper);
        pro_dis_vec(n2,1:m)=state_space(n2,1:m);
        pro_dis_vec(n2,m+1)=prob_visit;
    end
    
%     if ok==0
%         pro_dis_vec(n2,1:m)=S(n2,1:end);
%         pro_dis_vec(n2,m+1)=0;
%     end
    n2=n2+1;
end

%% One change for each input function

l=21;
DD=D(:,l-1);
DDP=D(:,l);


for i=1:m
    DDP(i)=DDP(i)+1;
    value_function(i)=ADP_doing_l(pro_dis_vec,r_hat,DD,DDP,C,N,lambda,mu,i,l);
    DDP(i)=DDP(i)-1;
end

%%

%% Evaluating associated value function to each input

c_value_function= value_function + mean_expected_revenue;

% for i=1:st_size(1)
%     if initial_state==S(i,:)
%        index_state=i;
%     end
% end
% phi_cal2(m,transition(n2,1:NT(1,2)))

ex_c_value_function=lambda'.*c_value_function+(1-lambda)'.*(phi_cal2(m,initial_state)'*r_hat);
re_c_value_function=phi_cal2(m,initial_state)'*r_hat+mean_expected_revenue;

C_value_function=[ex_c_value_function re_c_value_function];

[Reservation_revenue Decision_reservation]=max(C_value_function);


toc

