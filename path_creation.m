function [next_state] = path_creation(transition,combined_probability,f,a,s)

n1=size(transition);
n2=1;
n3=1;
%f=rand
while n2<=n1(1)
    combined(n2,1)=prod(combined_probability(n2,:));
    n2=n2+1;
end

%f=rand;
n3=1;
%transition

%% The selection of the transition by camparing the "f=rand" by the
%"combined_probability"
while n3 <= n1(1)
    y=sum(combined(1:n3));
    if f <= y
       next_state=transition(n3,1:end);
       break
    end
    n3=n3+1;
end

