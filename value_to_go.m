function [Next_Revenue] = value_to_go(transition,combined_probability,r_hat)

NT=size(transition);
n1=1;
n2=1;
n3=0;
sum_Next_Revenue=0;
%Discount_Factor=1;

while n1<= NT(1,1)
%    while n2<= NS
        
        %b=isequal(transition(n1,1:NT(1,2)),x(n2,1:NT(1,2)));
        
        y1=phi_cal2(NT(1,2),transition(n1,1:NT(1,2)))'*r_hat;
        y2=y1*prod(combined_probability(n1,1:NT(1,2)));
        sum_Next_Revenue=sum_Next_Revenue+y2;
%         if b==1
%             n3=n3+1;
%             n2;
%             betta=x(n2,end-1);
%             sum_Next_Revenue(1,n3)=x(n2,end-1)*prod(combined_probability(n1,1:NT(1,2)));
%         end
%        n2=n2+1;
%    end
    
    n1=n1+1;
    n2=1;
end
Next_Revenue=sum_Next_Revenue;