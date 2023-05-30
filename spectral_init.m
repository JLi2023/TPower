function [x1] = spectral_init(A,x,y_abs,s)
    
    [m,n] = size(A);
    % phase1: estimate support
    Y = create_Y(A,y_abs,1e8,0);
    diagY = diag(Y); [~,Omega] = maxk(abs(diagY),s); % maxY = max(diagY); j0 = diagY==maxY;

    v1 = svd_power(Y(Omega,Omega));
    xhat = zeros(n,1); xhat(Omega) = v1; 
    if norm(xhat+x)<norm(xhat-x), xhat = -xhat; end

%         % phase2: Truncated power method (optional)
%         ss = s;
%         [x0,~] = tpower(xhat,ss,Y); [x0,~] = proj_maxk(x0,s);
% %         [x0,~] = tpower_rand(xhat,ss,Y,A,y_abs,u,l); [x0,~] = proj_maxk(x0,s); % random sampling
%         if norm(x0+x)<norm(x0-x), x0 = -x0; end

    % phase3: refinement
    x1 = htp(A,xhat,x,y_abs,s);
    % fprintf('error %f\n', min(norm(w-x),norm(w+x))/norm(x));   
end

