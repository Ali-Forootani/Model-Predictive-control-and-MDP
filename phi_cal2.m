function [basis_val] = phi_cal2(m,s)

%[size_basis]= basis_cal(m);
n1=1;
n5=0;

    if n1==1
        phi(n1)=1;
    end
    
    for n1=2:(m+1)
        phi(n1)=s(n1-1);
    end 
   
%     for n2=1:m
%         for n3=n2:m
%             n5=n5+1;
%             r(n5)=s(n2)*s(n3); 
%         end
%      end
%     
%     phi(m+2:size_basis)=r;
    
    basis_val=phi';
  