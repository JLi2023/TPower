% power method
function max_sv = svd_power(M)

[m, n] = size(M);
x = randn(n,1);
for it = 1:100
   y = M*x;
   y = y/norm(y);
   if norm(x-y)/norm(x) < 1e-6
        break;
   end
   x = y;
end
   max_sv = x; 
end