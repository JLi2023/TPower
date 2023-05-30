function [x_rd,it] = htp(A,x0,x_true,y,s)

[m,n] = size(A);
mu = 0.95;
maxit = 20;
tol = 1e-4;
it=0;
stop=0;
%==========
%normAi = full(abs(sum(A'.*A',1)));
%cumul = disrand(normAi);
%========blockcircle
 
while ~stop
 it=it+1;
   err = min(norm(x0-x_true),norm(x0+x_true))/norm(x_true,'fro');
    %logerr=log(err);
    Ax = A*x0;
    y_sign=y.*sign(Ax);
    %%=======truncated HTP==========================
%     It_mag = abs(Ax);
%     T= find(It_mag > y/(1+gamma));
%     grad_i=A(T,:)'*(A(T,:)*x0-y_sign(T))/m;
    %%==========HTP==================================
    grad_i=A'*(A*x0-y_sign)/m;
    x1=x0-mu*grad_i;
    [~,num]=sort(abs(x1),'descend');
    S=num(1:s);
    As=A(:,S);
    Y=As'*As;
    z1=zeros(n,1);
    z1(S)=Y\(As'*y_sign);
%   z1=zeros(n,1);
%   z1(num(1:s))=x1(num(1:s)); 
    x0=z1;
%     fprintf('error %f\n', min(norm(x0-x_true),norm(x0+x_true)));
    %loghis(it)=logerr;
    %fprintf( 'Intermediate Soln:Error of recovered signal is %3.3e \n', err);
    if it>=maxit ||err<tol
         stop =1;
    end
   if it>=maxit
        stop =1;
   end
end
x_rd=x0;
end