function [b] = Pair_finding_ADP_dis(new_state_space,s)

NT=size(new_state_space);
n1=1;
n2=1;
n3=0;

while n1<= NT(1,1)
    %while n2<= NS
        
        x=isequal(new_state_space(n1,1:NT(1,2)),s);
        
        if x==0
            n3=n3+1;
            %b=0;
        %else
        %    b=1;
        end
        
%         if b==0
%             n3=n3+1;
%             new_state_visits(n3,1:m)=s;
%             n2;
%             %betta=x(n2,end-1);
%             %sum_Next_Revenue(1,n3)=x(n2,end-1)*prod(combined_probability(n1,1:NT(1,2)));
%         end
        %n2=n2+1;
    %end
    
    n1=n1+1;
    %n2=1;
end

if n3==NT(1,1)
    b=0;
else
    b=1;
end
%frequency_visit=n3;
%prob_visit=n3/num_exper;
