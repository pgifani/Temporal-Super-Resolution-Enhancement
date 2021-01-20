function x_BCS=bayesian_sparse_coding(y, A)
 

y =y(:,1);


N=size(A,2);


initsigma2 = std(y)^2/1e6;

[weights,used,~,~] = SBL(A,y,initsigma2,1e-5);

x_BCS = zeros(N,1);

x_BCS(used) = weights; 


