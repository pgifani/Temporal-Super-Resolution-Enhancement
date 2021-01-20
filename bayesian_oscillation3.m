function [ x_BCS,used,sigma2,errbars]=bayesian_oscillation3(y, A, select)

path(path, './Algorithm');

% x =x(:,select);
y =y(:,select);

% K=size(y,1);
% sigma=e;
N=size(A,2);
% % solve by BP
% x0 = A'*inv(A*A')*y;
% % take epsilon a little bigger than sigma*sqrt(K)
% epsilon =  sigma*sqrt(K)*sqrt(1 + 2*sqrt(2)/sqrt(K));                                                                                                              
% tic;
% x_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);
% t_BP = toc;
% fprintf(1,'BP number of nonzero weights: %d\n',sum(x_BP~=0));
% x_BP=x;


% solve by BCS
initsigma2 = std(y)^2/1e6;
% initsigma2 = 15;
% tic;
[weights,used,sigma2,errbars] = SBL(A,y,initsigma2,1e-5);
% t_BCS = toc;
% fprintf(1,'BCS number of nonzero weights: %d\n',length(used));
x_BCS = zeros(N,1);
% err = zeros(N,1);
x_BCS(used) = weights; 
% err(used) = errbars;

% reconstruction error
% E_BP = norm(x-x_BP)/norm(x);
% E_BCS = norm(x-x_BCS)/norm(x);

% Measure the error between the estimated parameter and the original one
% parameter_error=norm(x_BCS-x);
% dynamics_error=norm(A*x_BCS-y); 

% figure
% subplot(2,1,1); plot(x); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(a) Original Signal']);
% % subplot(3,1,2); plot(x_BP); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(b) Reconstruction with BP, K=' num2str(K)]);
% subplot(2,1,2); errorbar(x_BCS,err); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(c) Reconstruction with BCS, K=' num2str(K)]); box on;
% 
% % disp(['BP: ||I_hat-I||/||I|| = ' num2str(E_BP) ', time = ' num2str(t_BP) ' secs']);
% disp(['BCS: ||I_hat-I||/||I|| = ' num2str(E_BCS) ', time = ' num2str(t_BCS) ' secs']);