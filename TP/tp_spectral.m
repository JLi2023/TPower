function [x1] = tp_spectral(A,x,y_abs,s,u,l,b)
    
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
        [x0,~] = tpower_spectral(xhat,ss,Y,A,y_abs,1); [x0,~] = proj_maxk(x0,s);
%         [x0,~] = tpower_rand(xhat,ss,Y,A,y_abs,u,l); [x0,~] = proj_maxk(x0,s); % random sampling
        if norm(x0+x)<norm(x0-x), x0 = -x0; end

        % phase3: refinement
        w = htp(A,x0,x,y_abs,s);
        % fprintf('error %f\n', min(norm(w-x),norm(w+x))/norm(x));
        if norm(grad(A,y_abs,w)) < minnorm
            x1 = w;
            minnorm = norm(grad(A,y_abs,w));
            idx = l;
        end
    end
%     fprintf('Output the %d -th largest\n', idx);
end

function g = grad(A,y_abs,w)
g = A'*(A*w - y_abs.*sign(A*w));
end