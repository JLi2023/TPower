function [Y] = create_Y(A,y_abs,u,l)

[m,n] = size(A);
act_idx = find(y_abs<u & y_abs>l);
yy = zeros(m,1);
yy(act_idx) = y_abs(act_idx);
yy = diag(yy.^2);
Y = A'*yy*A;
% for i = 1:m
%     if y_abs(i) < u && y_abs(i) > l
%         Y = Y + y_abs(i)^2 * A(i,:)' * A(i,:);
%     end
% end
Y = Y/m;