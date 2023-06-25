function [combined_probability] = expectation_re(s,N,m,transition,mu,type)

%%
NT=size(transition);
n1=1;
n2=1;
n3=0;
n4=0;
n5=0;
n6=1;
n7=0;
n8=0;
n9=0;
n10=0;
n13=0;
n14=0;
n15=0;
n16=0;
NZ=0;
DZ=0;
BZ=0;
OR=0;
NZero=0;
NNega=0;
None=0;
transition;
probability=zeros(NT(1,1),m);

%%

while n6 <= m
    if s(n6)==0
        NZ=NZ+1;
        n7=n7+1;
        inNZ(1,n7)=n6;
    end
    if s(n6)~=0 
        n8=n8+1;
        DZ=DZ+1;
        inDZ(1,n8)=n6;
    end
    n6=n6+1;
end


while n1<= NT(1,1)
    val=transition(n1,1:end)-s;
    
     while n2 <= m 
        if val(n2)==0 & type~=6
            NZero=NZero+1;
            n3=n3+1;
            inNZero(1,n3)=n2;
        end        

        if val(n2)==-1
            n4=n4+1;
            NNega=NNega+1;
            inNNega(1,n4)=n2;
        end

         n2=n2+1;
     end
%     
    if type==1
        nothing(inDZ,1)=1-mu(inDZ,1);
    end
    
    if type==2
        nothing(inDZ,1)=1-mu(inDZ,1);
        nothing(inNZ,1)=1;
    end
    
    
    if type==4
        nothing(inNZ,1)=1;
    end

%      
     if NZero~=0 & type<=4
         probability(n1,[inNZero])=nothing([inNZero],1);
     end
%      
     if NNega~=0
         probability(n1,[inNNega])=mu([inNNega],1);
     end
%      
    n1=n1+1;
    n2=1;
    n3=0;
    n4=0;
    n5=0;
    None=0;
    NZero=0;
    NNega=0;
    
    clear inNone inNNega inNZero 
end
NNega;
%nothing;
 n13=0;
 probability;
 summation=0;
while n13 < NT(1,1)
    n13=n13+1;
    summation=prod(probability(n13,1:end))+summation;
end
summation;
combined_probability=probability;

