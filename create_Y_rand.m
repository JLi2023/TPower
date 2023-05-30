function [Y] = create_Y_rand(A,y_abs,p)

[m,n] = size(A);
act_idx = find(rand(m,1)<p);

yy = zeros(m,1);
yy(act_idx) = y_abs(act_idx);
yy = diag(yy.^2);
Y = A'*yy*A/m;

