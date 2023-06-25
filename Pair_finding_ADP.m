function [rep_ok] = Pair_finding_ADP(state_space,s)

NT=size(state_space);
n1=1;
n2=1;
n3=0;

while n1<= NT(1,1)
   b=isequal(state_space(n1,1:NT(1,2)),s);
   if b==1
      rep_ok=1;
      break
      %n3=n3+1;
      %n2;
   else
       rep_ok=0;
   end
   
    n1=n1+1;
end
%frequency_visit=n3;
%prob_visit=n3/num_exper;


