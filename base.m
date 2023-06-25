function [s] = base(m,N,n)
r=0;
qutient=1;
x=m;
while qutient~= 0
    qutient=floor(n/(N+1));
    remainder=mod(n,(N+1));
    s(1,x-r)=remainder;
    n=qutient;
    r=r+1;
end