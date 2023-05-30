function [x1,x10] = tpII(A,x,y_abs,s,u,l,b)
    
    [m,n] = size(A);
    % phase1: estimate support
    Y = create_Y(A,y_abs,u,l);
%     Y = create_Y_rand(A,y_abs,0.9);
    diagY = diag(Y); [~,maxb] = maxk(abs(diagY),b); % maxY = max(diagY); j0 = diagY==maxY;
    minnorm = 1e8;
    for l = 1:b
        j0 = maxb(l);
        ej0 = zeros(n,1); ej0(j0) = 1;
        [~,hatOmega] = maxk(abs(Y*ej0),s);
        v1 = svd_power(Y(hatOmega,hatOmega));
        xhat = zeros(n,1); xhat(hatOmega) = v1; 
        if norm(xhat+x)<norm(xhat-x), xhat = -xhat; end

        % phase2: Truncated power method (optional)
        ss = s;
        [x0,~] = tpower(xhat,ss,Y); [x0,~] = proj_maxk(x0,s);
%         [x0,~] = tpower_rand(xhat,ss,Y,A,y_abs,u,l); [x0,~] = proj_maxk(x0,s); % random sampling
        if norm(x0+x)<norm(x0-x), x0 = -x0; end
        
        if norm(grad(A,y_abs,x0)) < minnorm
            w = x0;
            minnorm = norm(grad(A,y_abs,x0));
            idx = l;
        end

        
        
    end
    % phase3: refinement
    x1 = htp(A,w,x,y_abs,s);
    % fprintf('error %f\n', min(norm(w-x),norm(w+x))/norm(x));
end

function g = grad(A,y_abs,w)
g = A'*(A*w - y_abs.*sign(A*w));
end