clear all
clc

%%

%load D_N6_m3.mat
%D=zeros(m,T);
%load Revenue_N6_m3.mat
%lambda=[0.4;0.4;0.3];
%mu=[0.2;0.3;0.5];

lambda_max=[0.55;0.5;0.3];
mu_min=[0.2;0.2;0.4];
lambda_min=[0.3;0.2;0.1];
mu_max=[0.5;0.6;0.6];

C=[0.9 1 1.1];
m=3;
T=40;
N=6;
%k=9;
%l=21;
D=zeros(m,T);
F=zeros(m,T);

n45=1;
%%
%D_h=D;


while n45 <= 20
    
    
%%

f_3=rand;

if f_3 <= 0.5
    D_h=D;
    F_h=F;
    h_1=randi([n45+2 n45+10],1,1);
    h_2=randi([h_1+2 h_1+15],1,1);
    l=h_2;
    ox_D=1;
    ox_F=0;
else
    F_h=F;
    D_h=D;
    h_1=randi([n45+2 n45+10],1,1);
    h_2=h_1;
    l=h_2;
    ox_D=0;
    ox_F=1;
end

%%

[Revenue_h S] = Exact_DP_booking_unbounded_bounded_function(N,m,T,D_h,F_h,lambda_max,lambda_min,mu_max,mu_min,C);

Revenue=Revenue_h;

%%

k=n45;

ok5=0;
st_size=size(S);
 while ok5==0
    h=randi([1,st_size(1)],1,1);
    initial_state=S(h,1:m);
    ok5=state_difference(initial_state,D(:,k));
 end

D_g=D(:,h_1:h_2);
F_g=F(:,h_1:h_2);
n53=0;
for i=h_1:h_2
    n53=n53+1;
    
    if ox_D==1 & ox_F==0
    
        if norm(D_g(:,n53)+F_g(:,n53),1)< N
            ox=1;
        else
            ox=0;
        end
    else
        if norm(D_g(:,n53)+F_g(:,n53),1)< N
            ox=1;
        else
            ox=0;
        end
    
    end
end


%%

st_size=size(S);
n=1;
n2=1;
n3=1;
num_exper=100;

%%

while l >= h_1 & ox==1

%% 
      
    
%% Main Simulation


while n <= num_exper
 
    state=initial_state;
    n1=k;
    Expected_Revenue=[];
    Expected_Revenue=0;
    path_state=[];
    
    while n1<=l-2
        f=rand;
        %state=state+D(:,n1)';
        D_k=D(:,n1);
        D_k_p=D(:,n1+1);
        
        F_k=F(:,k);
        F_k_p=F(:,k+1);
        
        n2=n2+1;
        %% Expected Revenue for the Duration of k-(l-2)
        state;
        Expected_Revenue=Expected_Revenue+Addition(state,C);
     
        %%
        
        if  norm(D_k+F_k,1)<N & norm(D_k_p+F_k_p,1)==N
        
                if sum(state) <= N
                    transition=D_k_p'+F_k_p';
                    state=transition;
                end
                
        end
                
    %%
    
    if norm(D_k+F_k,1)==N & norm(D_k_p+F_k_p,1)==N
            
            if state == D_k'+F_k'
                transition=D_k_p'+F_k_p';
                state=transition;
            end
            
    end
    
    %%
    
    if norm(D_k+F_k,1)==N & norm(D_k_p+F_k_p,1)< N
        
        for i=1:st_size(1)
            if state==S(i,:)
               index_state=i;
            end
        end
        n1;
        state;
        Decision=Revenue(index_state,m+2,n1);
        
        if Decision~=0 & Decision < m+1
                S_p=S(index_state,1:end)-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
        end
        
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
        
 
%%     
if norm(D_k_p+F_k_p,1) < N & norm(D_k+F_k,1) < N         
        
        state;
        if (sum(state)<N)
            
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
            n1;
            state;
            Decision=Revenue(index_state,m+2,n1);
            D_k';
            D_k_p';
            n1;
            if Decision~=0 & Decision < m+1
                S_p=S(index_state,1:end)-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
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
end         
       
        %%
        if norm(D_k+F_k,1)==0
            state;
            
            if (sum(state)<N)
            
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
            n1;
            state;
            Decision=Revenue(index_state,m+2,n1);
            D_k';
            D_k_p';
            n1;
            if Decision~=0 & Decision < m+1
                S_p=S(index_state,1:end)-D_k';
                [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                state_path=path_creation(transition,combined_probability,f);
                state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
            end
            
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
           
        end
        %%
        
        if (sum(state)==N)
            
            D_k;
            D_k_p;
            n1;
            state;
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
            if norm(D_k+F_k,1) == 0
                Decision=Revenue(index_state,m+2,n1);
                
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
               n1;
               index_state;
               Decision=Revenue(index_state,m+2,n1);
               index_state;
                if Decision~=0 & Decision < (m+1)
                    S_p=S(index_state,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,Decision,lambda,mu,D_k,D_k_p);
                    state_path=path_creation(transition,combined_probability,f);
                    state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
                end
                if Decision == (m+1)
                    S_p=S(index_state,1:end)-D_k';
                    [lambda mu] = lambda_mu_calculation(S_p,N,m,lambda_max,lambda_min,mu_max,mu_min);
                    [transition_ca combined_probability_ca]=stateanalysis_re_bo_ver_ca_3(S_p,N,m,mu,D_k,D_k_p);
                    state_path=path_creation_re(transition_ca,combined_probability_ca,f);
                    state=shifted_transition_unbounded(state_path,D_k_p,F_k_p,m);
                end
                
                if Decision == m+2
                ok_2=0;
                er=0;
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
        end
        state;
        n1=n1+1;
        path_state(n1,1:m)=state;
       
    end
    state_visits(n,1:m)=state;
    Expected_Revenue_Experiment(n)=Expected_Revenue;
    %n2=0;
    %"====================="
    n=n+1;
    %path_state=[];
end

%% mean expected revenue for the time duration of k-(l-2)
mean_expected_revenue=sum(Expected_Revenue_Experiment)/num_exper;

%% Computing the probability distribution for the time slot "l-1"

n2=1;
%num_exper=200;
while n2 <= st_size(1)
    D_k = D(:,l-1);
    F_k = F(:,l-1);
    [ok] = state_difference_unbounded(S(n2,1:end),D_k,F_k);
    if ok==1
        [frequency_visit prob_visit]=Pair_finding_dis(state_visits,S(n2,1:end),num_exper);
        pro_dis_vec(n2,1:m)=S(n2,1:end);
        pro_dis_vec(n2,m+1)=prob_visit;
    end
    
    if ok==0
        pro_dis_vec(n2,1:m)=S(n2,1:end);
        pro_dis_vec(n2,m+1)=0;
    end
    n2=n2+1;
end


%% One change for each input function

DD=D(:,l-1);
DDP=D(:,l);
FF=F(:,l-1);
FFP=F(:,l);

if ox_D==1 & ox_F==0
for i=1:m
    DDP(i)=DDP(i)+1;
    value_function(i) = DP_doing_l_online_unbounded(pro_dis_vec,Revenue,DD,DDP,FF,FFP,C,N,lambda_max,lambda_min,mu_max,mu_min,i,l);
    DDP(i)=DDP(i)-1;
end
end

if ox_D==0 & ox_F==1
for i=1:m
    FFP(i)=FFP(i)+1;
    value_function(i) = DP_doing_l_online_unbounded(pro_dis_vec,Revenue,DD,DDP,FF,FFP,C,N,lambda_max,lambda_min,mu_max,mu_min,i,l);
    FFP(i)=FFP(i)-1;
end
end


%% Evaluating associated value function to each input
c_value_function = value_function + mean_expected_revenue;

lambda_acc=[0.6;0.58;0.55];
mu_acc=[0.2;0.2;0.3];

%ex_c_value_function=lambda_acc'.*c_value_function+(1-lambda_acc)'.*Revenue_h(index_state,m+1,k);
%re_c_value_function=Revenue(index_state,m+1,k);
C_value_function = [c_value_function]

%C_value_function = [ex_c_value_function]

[Reservation_revenue Decision_reservation] = max(C_value_function);



%% Accepting or Rejecting from customer side

%D_h=D;
%if Decision_reservation < m+1
    %if f1 <= lambda(Decision_reservation)
    if ox_D==1 & ox_F==0
    D(Decision_reservation,l)=D(Decision_reservation,l)+1;
    %end
    [Revenue S] = Exact_DP_booking_unbounded_bounded_function(N,m,T,D,F,lambda_max,lambda_min,mu_max,mu_min,C);
    l=l-1;
    end
    
    if ox_F==1 & ox_D==0
       F(Decision_reservation,l)=F(Decision_reservation,l)+1;
       [Revenue S] = Exact_DP_booking_unbounded_bounded_function(N,m,T,D,F,lambda_max,lambda_min,mu_max,mu_min,C);
        l=l-1;
    end
    
%else
%    D=D_h;
    %D(:,l)=D(:,l);
    %[Revenue S] = Exact_DP_Booking_10(N,m,T,D,lambda_max,lambda_min,mu_max,mu_min,C);
%end

end

%%

if ox==1
 
for i=1:st_size(1)
    if initial_state == S(i,:)
       index_state=i;
    end
end

ex_c_value_function=lambda_acc'.*c_value_function+(1-lambda_acc)'.*Revenue_h(index_state,m+1,k);

C_value_function_2 = [ex_c_value_function Revenue_h(index_state,m+1,k)];

[Reservation_revenue Decision_reservation] = max(C_value_function_2);

if Decision_reservation ~= m+1
    
    D_h=D;
    F_h=F;
    
    %if ox_D==1 & ox_F==0
    %    D_h=D;
    %    F_h=F;
    %end
    
    %if ox_F==1 & ox_D==0
    %    D_
    %end
    
    %if f1 <= lambda(Decision_reservation)
    %D(Decision_reservation,l)=D(Decision_reservation,l)+1;
    %end
%    [Revenue S] = Exact_DP_Booking_10(N,m,T,D,lambda_max,lambda_min,mu_max,mu_min,C);
%    l=l-1;
%else
%    D=D_h;
    %D(:,l)=D(:,l);
    %[Revenue S] = Exact_DP_Booking_10(N,m,T,D,lambda_max,lambda_min,mu_max,mu_min,C);
end


end

reservation_request(n45,1:2)=[h_1 h_2];

n45=n45+1;

end


%% 

[Revenue_h S] = Exact_DP_booking_unbounded_bounded_function(N,m,T,D_h,F_h,lambda_max,lambda_min,mu_max,mu_min,C);

%%

