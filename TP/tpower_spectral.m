function [w,Omega] = tpower_spectral(z,ss,Y,A,y_abs,random_sample)

[m,n] = size(A);
maxit = 10;
p = 0.8;
for i = 1:maxit
    if random_sample 
        SS = find(rand(m,1)<p);
        Y = create_Y(A(SS,:),y_abs(SS),1e8,0);
    end
    
    zprev = z;
    ww = Y*z;
    [~,Omega] = maxk(abs(ww),ss);
    v1 = svd_power(Y(Omega,Omega));
    z = zeros(n,1); z(Omega) = v1; 
    if norm(z-zprev)<1e-4
        break
    end
end
w = z;
