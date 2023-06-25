function [value_function_ADP] = ADP_doing_l(pro_dis_vec,r_hat,DD,DDP,C,N,lambda,mu,a,l)
state_space=pro_dis_vec(:,1:end-1);
NT=size(state_space);
m=NT(1,2);
n1=1;
n2=1;
n3=0;
D_k=DD;
D_k_p=DDP;
value_function_ADP=0;

%x=Revenue(1:end,1:end,l);
while n1<= NT(1,1)
    %while n2<= NS
    if pro_dis_vec(n1,m+1)~=0
        state=pro_dis_vec(n1,1:m);
        S_p=state-D_k';
        [transition combined_probability]=stateanalysis_bo_ver_3(S_p,N,m,a,lambda,mu,D_k,D_k_p);
        transition_2 = shifted_transition(transition,D_k_p,m);
        Rev = ADP_r_evaluation(transition_2,r_hat,pro_dis_vec,NT(1,1),combined_probability);
        %Rev=Pair_finding(transition_2,x,NT(1,1),combined_probability);
        value_function_ADP=pro_dis_vec(n1,m+1)*(Addition(state,C)+Rev)+value_function_ADP;
    end
        %b=isequal(state_visits(n1,1:NT(1,2)),s);
        
        %if b==1
        %    n3=n3+1;
        %    n2;
            %betta=x(n2,end-1);
            %sum_Next_Revenue(1,n3)=x(n2,end-1)*prod(combined_probability(n1,1:NT(1,2)));
        %end
        %n2=n2+1;
    %end
    
    n1=n1+1;
    %n2=1;
end
%frequency_visit=n3;
%value_function=n3/num_exper;
