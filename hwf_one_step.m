function [x1,x10] = hwf_one_step(A,x,y_abs,s,b)
    [m,n] = size(A);
    minnorm = 1e80;
    Y = create_Y(A,y_abs,1e8,0);
    diagY = diag(Y); [~,maxb] = maxk(abs(diagY),b); % maxY = max(diagY); j0 = diagY==maxY;
    for l = 1:b
        j0 = maxb(l);
        ej0 = zeros(n,1); ej0(j0) = 1;
        [~,hatOmega] = maxk(abs(Y*ej0),s);
        v1 = svd_power(Y(hatOmega,hatOmega));
        xhat = zeros(n,1); xhat(hatOmega) = v1; 
        if norm(xhat+x)<norm(xhat-x), xhat = -xhat; end
        w = htp(A,xhat,x,y_abs,s); 
%         [w] = sparse_stochasticADM(xhat,x,s,A,y_abs,0.6,100,1e-5);
%         %% start descent
%         maxit = 50;
%         mu = 1;
%         gamma = 0.7;
%         tol2 = 1e-4;
%         w = xhat;
%         for t = 1:maxit
%             It_act = A*w;
%             It_mag = abs(It_act);
%             It = find(It_mag > y_abs/(1+gamma));
%             sum_TAF = zeros(n,1);
%             for i=1:length(It)
%                 ii = It(i); 
%                 sum_TAF = sum_TAF + (It_act(ii,1) - y_abs(ii)*(It_act(ii,1))/It_mag(ii,1))*A(ii,:)'; %(n x 1)
%             end
%             grad_x = (mu/m)*sum_TAF;
%             arg_TAF = w - grad_x;
%             w = truncated_AF(arg_TAF,s);
%             err_hist(t+1,1) = norm(y_abs-abs(A*w))/norm(y_abs);
%             err_hist(t+1,2) = min(norm(w-x),norm(w+x))/norm(x);
%             if err_hist(t+1,2)<tol2
%                 break;
%             end
%         end
        
        if norm(grad(A,y_abs,w)) < minnorm
            x1 = w;
            minnorm = norm(grad(A,y_abs,w));
            idx = l;
        end
        
        if l == 1
            x10 = w;
        end
    end
%     fprintf('Output the %d -th largest\n', idx);
end


function g = grad(A,y_abs,w)
    [m,n] = size(A);
    Aw = A*w;
    b = (Aw.^2 - y_abs.^2).*Aw;
    g = A'*b/m;
end
