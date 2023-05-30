function [x1] = tp_sam(A,x,y_abs,s,u,l,b)
    
    [m,n] = size(A);
    % phase1: estimate support
    Y = create_Y(A,y_abs,u,l);
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
        if norm(x0+x)<norm(x0-x), x0 = -x0; end

        % phase3: refinement
        [w] = sparse_stochasticADM(x0,x,s,A,y_abs,0.6,100,1e-5);
        if norm(grad(A,y_abs,w)) < minnorm
            x1 = w;
            minnorm = norm(grad(A,y_abs,w));
        end
    end
end

function g = grad(A,y_abs,w)
g = A'*(A*w - y_abs.*sign(A*w));
end