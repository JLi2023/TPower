function [w,Omega] = tpower_rand(z,ss,Y,A,y_abs,u,l)

[m,n] = size(A);
random_sample = 1;
p = 0.9;
n = length(z);
maxit = 50;
for i = 1:maxit
    if random_sample 
        SS = find(rand(m,1)<p);
        Y = create_Y(A(SS,:),y_abs(SS),u,l);
    end
    zprev = z;
    ww = Y*z;
    [~,Omega] = maxk(abs(ww),ss);
    z = zeros(n,1); z(Omega) = ww(Omega); z = z/norm(z);
    if norm(z-zprev)<1e-4
        break
    end
end
w = z;