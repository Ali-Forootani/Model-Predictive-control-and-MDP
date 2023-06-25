function [s] = binbase(m,n)
r=0;
qutient=1;
x=m;
while qutient~= 0
    qutient=floor(n/(2));
    remainder=mod(n,(2));
    s(1,x-r)=remainder;
    n=qutient;
    r=r+1;
end