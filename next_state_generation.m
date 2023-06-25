function [next_state]= next_state_generation(m,N,Revenue,lambda_max,lambda_min...
,mu_max,mu_min,C,D_k,D_k_p,F_k,F_k_p,T,state)
        
        %%
        
        f=rand;
        %state=state+D(:,n1)';
        %D_k=D(:,n1);
        %D_k_p=D(:,n1+1);
        
        %F_k=F(:,k);
        %F_k_p=F(:,k+1);
        
        S = Revenue(:,1:m);
        st_size = size(S);
        
        less_less=0;
        equal_less=0;
        curr_demand=0;
        full_system=0;
        
        %n2=n2+1;
        %% Expected Revenue for the Duration of k-(l-2)
        %state;
        %Expected_Revenue=Expected_Revenue+Addition(state,C);
     
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
    
    if norm(D_k+F_k,1)==N & norm(D_k_p+F_k_p,1)< N & equal_less==0
        
        equal_less=1;
        
        
        for i=1:st_size(1)
            if state==S(i,:)
               index_state=i;
            end
        end
        state;
        Decision=Revenue(index_state,m+2);
        
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
                S_p=S(index_state,1:end)-D_k';
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
        
 
%%     
if norm(D_k_p+F_k_p,1) < N & norm(D_k+F_k,1) < N & less_less==0 & equal_less==0
        
        less_less=1;
        
        state;
        if (sum(state)<N)
            
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
            
            state;
            Decision=Revenue(index_state,m+2);
            D_k';
            D_k_p';
            
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
end         
       
        %%
        if norm(D_k+F_k,1)==0 & less_less==0 & equal_less==0 & curr_demand==0
            
            curr_demand=1;
            
            state;
            
            if (sum(state)<N)
            
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
           
            state;
            Decision=Revenue(index_state,m+2);
            D_k';
            D_k_p';
            
            
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
           
        end
        %%
        
        if (sum(state)==N) & less_less==0 & equal_less==0 & curr_demand==0
            
            full_system=1;
            
            D_k;
            D_k_p;
            
            state;
            for i=1:st_size(1)
                if state==S(i,:)
                    index_state=i;
                end
            end
            if norm(D_k+F_k,1) == 0
                
                Decision=Revenue(index_state,m+2);
                
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
               
               index_state;
               Decision=Revenue(index_state,m+2);
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
                S_p=S(index_state,1:end)-D_k';
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
        end
        
        next_state=state;
        %n1=n1+1;
        %path_state(n1,1:m)=state;
       