function [transition combined_probability] = stateanalysis_re_bo_ver_ca_3(s,N,m,mu,D_k,D_k_p)
%% This function gets a state "s", "N", "m", "u_i=a", "lambda's", "mu's"
% Inputs: "s", "N", "m", "u_i=a", "lambda's", "mu's"
%
% Outputs: matrix "transition" and the matrix "combined_probability" of a going from
% the state "x(k)" to satate "x(k+1)", each row of the transition matrix
% has a probability!
%
%In general we have "11" types of the states, e.g. state s=[1 1] has
% different number of transitions campare to state s=[0 0] so in some way
% we have to find out how many transitions exist for a given state**
% we have to know that we are on the boundry or not**
% "binbase": It is a written function which help us to find out "x(k)"
% from a "x(k-1)" if action has been taken into account, e.g. x(k-1)=[1 1]
% x(k) could be [1 1], [0 0], [1 0], [0 1], [2 0] under the action
% "u_i=c_1"
% "expectation": It is a written function which gives the probability of a
% transtion from a state "x(k)" to the state "x(k+1)".  

%=====================================
%% auxilary variables
n2=1;
NZ=0;
n3=0;
n5=0;
n6=0;
n4=0;
n7=0;
n8=0;
n9=0;
n10=0;
n11=0;
n12=0;
n13=0;
n14=0;
n15=1;
n16=0;
n17=1;
n18=0;
n19=1;
n20=0;
n21=1;
n22=0;
n23=1;
n24=0;
n25=1;
n27=1;
n29=1;
n30=1;
BZ=0;
DZ=0;
OR=0;
clear transition tr
%% Suppose that state s=[1 0 0 3 4]
% "NZ": number of zero elements in a vector "s", for the example above NZ=2
% and there is no customer corresponding to price "c_2" and "c_3".

% "DZ": number of non-zero elements in a vector "s", e.g. for the above
% example BZ=3 and there is 1 customer in price "c_1", 3 customers in price
% "c_4", and 4 customers in price "c_5"

% "BZ": number of non-zero elements in a vector "s", and non-zero elements
% are not in the group related to action "u_i=a"

% "OR": number of zero elements in a vector "s", and zero elements are not
% in the group related to action "u_i=a"

while n2 <= m
    if s(n2)==0
        NZ=NZ+1;
        n3=n3+1;
        inNZ(1,n3)=n2;
    end
    if s(n2)~=0 
        n4=n4+1;
        DZ=DZ+1;
        inDZ(1,n4)=n2;
    end
%     if s(n2)~=0 & n2~= a
%         n9=n9+1;
%         BZ=BZ+1;
%         inBZ(1,n9)=n2;
%     end
%     if s(n2)==0 & n2~= a
%         n24=n24+1;
%         OR=OR+1;
%         inOR(1,n24)=n2;
%     end
    n2=n2+1;
end

% Having to know how many non-zero elements, zero elements, zero customer
% corresponding to the "u_i=a" help us to figure out for a given state how
% many transition could be possible
%%

%In general we are on the "boundry" or not!
% so at first we check "sum(s)<N" 

% "type=1": state like s=[1 2 3], N=7, m=3, and action "u_i=c_2"
% e.g. here we have 3*2^(2) transitions

%===============

% "type=2": state like s=[1 2 0 2], N=7, m=4, and action "u_i=c_2"
% e.g. here we have 3*2^(2)=12 transitions, since "x_3=0"!
%===============

% "type=3": state like s=[3 0 0 0], N=7, m=4, action "u_i=c_1"
% e.g. here we have 3 transitions, since "x_2=0,x_3=0,x_4=0"!

%===============
% "type=4": state like s=[0 0 0 0], N=7, m=4, action "u_i=c_1"
% e.g. here we have 2 transitions, since we are at the origin

%===============
% "type=5": state like s=[0 2 1 1], N=7, m=4, action "u_i=c_1", "x_1=0"!!!!
% e.g. here we have 2*2^(4-1)=16, since "x_1=0" and we choose "u_i=c_1"

if sum(s)<=N
    if NZ==0
        transition=zeros(2^(m),m);
        type=1;
        while n5 <= 2^m-1
            tr=-1*binbase(m,n5);
            n5=n5+1;
            transition(n5,1:end)=s+tr;
        end        
    end
    if NZ~=0 & NZ<=(m-1)
            transition=zeros(2^(m-NZ),m);
            type=2;
            while n7 <= 2^(m-NZ)-1
                tr=-1*binbase(m-NZ,n7);
                n7=n7+1;
                transition(n7,inDZ)=s(1,inDZ)+tr;
            end
    end
%     if NZ~=0 & NZ==(m-1)
%               transition=zeros(2^(m-NZ),m);
%               type=3;
%               while n10 <= 2^(m-NZ)-1
%                 tr=-1*binbase(m-NZ,n10);
%                 n10=n10+1;
%                 transition(n10,inDZ)=s(1,inDZ)+tr;
%               end
%     end
    if s==0
                type=4;
                transition=zeros(1,m);
                transition(1,1:end)=s;
                %transition(2,a)=s(a)+1;
    end
    
%     if s(a)==0 & BZ~=0
%         transition=zeros(2*2^(BZ),m);
%         type=5;
%         while n12 <= 2^(BZ)-1
%             tr=-1*binbase(BZ,n12);
%             n12=n12+1;
%             transition(n12,inDZ)=s(1,inDZ)+tr;
%         end
%         while n13 <= 2^(BZ)-1
%             tr=-1*binbase(BZ,n13);  
%             n13=n13+1;
%             n12=n12+1;
%             transition(n12,inDZ)=s(1,inDZ)+tr;
%             transition(n12,a)=s(a)+1;
%         end
%     end
    type;
    combined_probability=expectation_re(s,N,m,transition,mu,type);
end



%%
% 
% n25=1;
% n28=1;
% 
% tr_size=size(transition);
%     while n25 <= tr_size(1)
%         sub_tract=transition(n25,:)'-D_k_p;
%         while n28 <= m 
%             if transition(n25,n28)< D_k_p(n28)
%                 transition(n25,n28)=D_k_p(n28);
%             end
%             n28=n28+1;
%         end
%         n25=n25+1;
%         n28=1;
%     end




