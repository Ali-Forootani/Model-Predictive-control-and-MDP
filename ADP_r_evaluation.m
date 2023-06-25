function [Next_Revenue] = ADP_r_evaluation(transition,r_hat,pro_dis_vec,NS,combined_probability)

NT=size(transition);
n1=1;
n2=1;
n3=0;
m=NT(1,2);
%while n1<= NT(1,1)
    while n2<= NT(1,1)
        
        %b=isequal(transition(n1,1:NT(1,2)),pro_dis_vec(n2,1:NT(1,2)));
        
        %if b==1
            n3=n3+1;
            n2;
            %betta=pro_dis_vec(n2,end);
            basis_func=phi_cal2(m,transition(n2,1:NT(1,2)));
            x_valu=basis_func'*r_hat;
            sum_Next_Revenue(1,n3)=x_valu*prod(combined_probability(n1,1:NT(1,2)));
        %end
        n2=n2+1;
    end
    
    
%    n1=n1+1;
%    n2=1;
%end
Next_Revenue=sum(sum_Next_Revenue);