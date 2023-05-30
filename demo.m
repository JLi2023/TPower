clc;clear;
addpath(genpath(pwd))
%% signal model
n = 1000;
s_set = [20]; slen = length(s_set);
m_set = 100:100:1000; mlen = length(m_set);
result = zeros(slen,mlen,6);
% other parameters
u = 5;
l = 0.8;
b = 20; % multistart: for TP and HWF

repeat = 50;
act = [1,0,0,0,0,0];

%% TP
if act(1)
    fprintf('Running TP......\n')
    rng('default');
    for i = 1:slen
        s = s_set(i);
        % generate signal
        [x,xsupp] = gen_signal(n,s);
        for j = 1:mlen
            count = 0;
            m = m_set(j);
            for ii = 1:repeat
                fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
                % generate measurements
                [y_abs,y_ph,A] = measure_signal(m,x);
                
                [w] = tp(A,x,y_abs,s,u,l,b);

                if min(norm(w-x),norm(w+x))/norm(x) < 1e-2
                    count = count + 1;
                end
            end
            result(i,j,1) = count;
        end
    end
end


% %% Hadamard
% if act(2)
%     fprintf('Running Hadamard......\n')
%     rng('default');
%     for i = 1:slen
%         s = s_set(i);
%         % generate signal
%         [x,xsupp] = gen_signal(n,s);
%         for j = 1:mlen
%             count = 0;
%             m = m_set(j);
%             for ii = 1:repeat
%                 fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
%                 % generate measurements
%                 [y_abs,y_ph,A] = measure_signal(m,x);
% 
%                 [w] = hwf_one_step(A,x,y_abs,s,b);
% %                 fprintf('error %f\n', min(norm(w-x),norm(w+x))/norm(x));
% 
%                 if min(norm(w-x),norm(w+x))/norm(x) < 1e-2
%                     count = count + 1;
%                 end
%             end
%             result(i,j,2) = count;
%         end
%     end
% end
% 
% %% tp + spectral
% if act(3)
%     fprintf('Running TP + spectral......\n')
%     rng('default');
%     for i = 1:slen
%         s = s_set(i);
%         % generate signal
%         [x,xsupp] = gen_signal(n,s);
%         for j = 1:mlen
%             count = 0;
%             m = m_set(j);
%             for ii = 1:repeat
%                 fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
%                 % generate measurements
%                 [y_abs,y_ph,A] = measure_signal(m,x);
% 
%                 %% initialization
%                 % phase1: estimate support
%                 Y = create_Y(A,y_abs,u,l);
%                 diagY = diag(Y); maxY = max(diagY); j0 = find(diagY==maxY);
%                 ej0 = zeros(n,1); ej0(j0) = 1;
%                 [~,hatOmega] = maxk(abs(Y*ej0),s);
%                 v1 = svd_power(Y(hatOmega,hatOmega));
%                 xhat = zeros(n,1); xhat(hatOmega) = v1; 
%                 if norm(xhat+x)<norm(xhat-x), xhat = -xhat; end
% 
%                 % phase2: Truncated power method (optional)
%                 ss = s;
%                 [x0,~] = tpower_spectral(xhat,ss,Y); [x0,hatOmega] = proj_maxk(x0,s);
%                 if norm(x0+x)<norm(x0-x), x0 = -x0; end
% 
%                 % phase3: refinement
%                 w = htp(A,x0,x,y_abs,s);
%                 fprintf('error %f\n', min(norm(w-x),norm(w+x))/norm(x));
% 
%                 if min(norm(w-x),norm(w+x))/norm(x) < 1e-2
%                     count = count + 1;
%                 end
%             end
%             result(i,j,3) = count;
%         end
%     end
% end
% 
% %% CoPRAM
% if act(4)
%     rng('default');
%     fprintf('Running CoPRAM......\n')
%     for i = 1:slen
%         s = s_set(i);
%         % generate signal
%         [x,xsupp] = gen_signal(n,s);
%         for j = 1:mlen
%             count = 0;
%             m = m_set(j);
%             for ii = 1:repeat
%                 fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
%                 % generate measurements
%                 [y_abs,y_ph,A] = measure_signal(m,x);
% 
%                 [x1,err_hist1,C1,x1_init] = CoPRAM(y_abs,A,s,20,1e-3,1e-4,x);
%                 if min(norm(x1-x),norm(x1+x))/norm(x) < 1e-2
%                     count = count + 1;
%                 end
%             end
%             result(i,j,4) = count;
%         end
%     end
% end
% 
% 
% %% ThWF
% if act(5)
%     rng('default');
%     fprintf('Running ThWF......\n')
%     for i = 1:slen
%         s = s_set(i);
%         % generate signal
%         [x,xsupp] = gen_signal(n,s);
%         for j = 1:mlen
%             count = 0;
%             m = m_set(j);
%             for ii = 1:repeat
%                 fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
%                 % generate measurements
%                 [y_abs,y_ph,A] = measure_signal(m,x);
%                 
%                 y_twf = y_abs.^2;  
%                 [x3] = Thresholded_WF(y_twf,A,s,200,0,0,x);
%                 if min(norm(x3-x),norm(x3+x))/norm(x) < 1e-2
%                     count = count + 1;
%                 end
%             end
%             result(i,j,5) = count;
%         end
%     end
% end
% 
% %% SPARTA
% if act(6)
%     rng('default');
%     fprintf('Running SPARTA......\n')
%     for i = 1:slen
%         s = s_set(i);
%         % generate signal
%         [x,xsupp] = gen_signal(n,s);
%         for j = 1:mlen
%             count = 0;
%             m = m_set(j);
%             for ii = 1:repeat
%                 fprintf('s:%d\t m:%d\t repeat:%d\n',s,m,ii);
%                 % generate measurements
%                 [y_abs,y_ph,A] = measure_signal(m,x);
%                 
%                 y_twf = y_abs.^2;  
%                 [x4] = SparTAF(y_abs,A,s,50,1e-3,1e-4,x);
%  
%                 if min(norm(x4-x),norm(x4+x))/norm(x) < 1e-2
%                     count = count + 1;
%                 end
%             end
%             result(i,j,6) = count;
%         end
%     end
% end
% 
% 
% 
% result = result/repeat;
% fprintf('================END================\n');
% 



%% plot
% result1 = result(:,:,1);
% result2 = result(:,:,2);
% % result3 = result(:,:,3);
% % result4 = result(:,:,4);
% % result5 = result(:,:,5);
% % result6 = result(:,:,6);
% figure('DefaultAxesFontSize',15)
% hold on
% grid minor
% plot(result1,'r-o');
% plot(result2,'b-*');
% % plot(result3,'m-*');
% % plot(result4,'g-.');
% % plot(result5,'y-.');
% % plot(result6,'k-.');
% 
% % legend('TP','Hadamard','TP+spectral','CoPRAM','ThWF','SPARTA','Location', 'Best')
% legend('TP-MR','HWF-MR','Location', 'Best')
% 
% 
% title('Fixed sparsity s = 20, n = 1000, repeat = 50')
% ylabel('Success rate');
% xticks(1:10);
% xticklabels(100:100:1000);
% xlabel('Sampling number m');



