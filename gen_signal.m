function [x,supp] = gen_signal(n,s)

x = zeros(n,1);
l = rand(n,1);
[~,supp] = maxk(l,s);
x(supp) = randn(s,1);
% normalization for simplicity
x = x/norm(x);
