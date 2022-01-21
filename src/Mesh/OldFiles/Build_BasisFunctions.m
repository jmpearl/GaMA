clear all 
close all
clc

syms u v w

phi = {u,v,w,1-u-w-v};
BF = 1;
d = [2,2,1,1];
n=sum(d);

for j = 1:4
for i = 1:d(j)
    
    BF = BF * (n*phi{j}-i+1)/i;
    
end
end


int(int(int(BF,w,[0,1-u-v]),v,[0,1-u]),u,[0,1])


% ----6th-----
%  -4/1200 +...
%  12*(1/350) + 12*(-1/280) + 6*(1/210) +...
%  4*3*(0) + 4*6*(1/280) + 4*(-3/560) +...
%  4*(3/140) + 6*(0); 