function [frequency_visit prob_visit] = Pair_finding_dis(state_visits,s,num_exper)

NT=size(state_visits);
n1=1;
n2=1;
n3=0;
while n1<= NT(1,1)
    %while n2<= NS
        
        b=isequal(state_visits(n1,1:NT(1,2)),s);
        
        if b==1
            n3=n3+1;
            n2;
            %betta=x(n2,end-1);
            %sum_Next_Revenue(1,n3)=x(n2,end-1)*prod(combined_probability(n1,1:NT(1,2)));
        end
        %n2=n2+1;
    %end
    
    n1=n1+1;
    %n2=1;
end
frequency_visit=n3;
prob_visit=n3/num_exper;
