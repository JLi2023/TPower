function [w] = hwf(A,y_abs,x)

    [m,n] = size(A);
    alpha = 0.001;
    eta = 0.1;
    maxit = 100;
    theta = sqrt(sum(y_abs.^2)/m);

    Y = A'*diag(y_abs.^2)*A;
    diagY = diag(Y); maxY = max(diagY); j0 = find(diagY==maxY);
    % initialization
    u = alpha * ones(n,1);
    v = alpha * ones(n,1);
    u(j0) = sqrt(theta/sqrt(3) + alpha^2);
    
    
    % iteration
    for i = 1:maxit
        w = u.*u-v.*v;
        fprintf('error: %f\n',min(norm(w-x),norm(w+x)));
        u = u.*(ones(n,1) - 2*eta*grad(A,y_abs,w));
        v = v.*(ones(n,1) + 2*eta*grad(A,y_abs,w));
    end
    
end 
    
    

function g = grad(A,y_abs,w)
    [m,n] = size(A);
    Aw = A*w;
    b = (Aw.^2 - y_abs.^2).*Aw;
    g = A'*b/m;
end