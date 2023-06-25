function [sample_space]= state_space_generation(m,N,sample_length,D,T)

Time_step=1;
n62=1;
%sample_length=1000;


while Time_step<=T
while n62<= sample_length
    g=(N-1)*rand(m,1);
    g=round(g);
    [ok] = state_difference(g',D(1:m,Time_step));
    if sum(g)<=N & ok==1
        sample_space(n62,1:m,Time_step)=g';
        n62=n62+1;
    end
end
Time_step=Time_step+1;
n62=1;
end