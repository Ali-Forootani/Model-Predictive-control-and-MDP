function [ok] = state_difference_unbounded(s,D_k,F_k)

sub_tract=s'-D_k-F_k;
n=0;
sub_size=size(sub_tract);

for i=1:sub_size(1)
    if sub_tract(i)>= 0
    n=n+1;    
    end
end

if n==sub_size(1)
    ok=1;
else
    ok=0;
end