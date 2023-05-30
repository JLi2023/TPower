function [w,Omega] = tpower(z,ss,Y)

n = length(z);
maxit = 50;
for i = 1:maxit
    zprev = z;
    ww = Y*z;
    [~,Omega] = maxk(abs(ww),ss);
    z = zeros(n,1); z(Omega) = ww(Omega); z = z/norm(z);
    if norm(z-zprev)<1e-4
        break
    end
end
w = z;