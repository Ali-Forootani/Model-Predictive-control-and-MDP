function [transition_2] = shifted_transition(transition,D_k_p,m)

tr_size=size(transition);
for n1=1:tr_size(1)
    y(n1,1:m)=D_k_p';
end
transition_2=transition+y;
% u=transition+y;
% n1=1;
% 
% while n1 <= tr_size(1)
%     if norm(u(n1,:),1)>N
%         state_status=0;
%         break
%     else
%         state_status=1;
%         n1=n1+1;
%     end
% end