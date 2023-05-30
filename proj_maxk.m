function [y,S] = proj_maxk(x,k)

y = zeros(length(x),1);
[~,S] = maxk(abs(x),k);
y(S) = x(S);
