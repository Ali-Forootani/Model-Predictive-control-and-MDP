function [next_state] = path_creation_re(transition_re,combined_probability_re,f)

n1=size(transition_re);
n2=1;
n3=1;
%f=rand
while n2<=n1(1)
    combined(n2,1)=prod(combined_probability_re(n2,:));
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
       next_state=transition_re(n3,1:end);
       break
    end
    n3=n3+1;
end

