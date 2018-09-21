function [r S] = monteCarlo(S0,mu,Q,T,N,nPaths)
% r = nPaths x n, mean return for each scenario(path)
rho=corrcov(Q);
sigmas=sqrt(diag(Q));
n=size(mu,1);
dt=T/N;
L=chol(rho,'lower');
S=zeros(N+1,n,nPaths);
S(1,:,:)=S0'.*ones(n,nPaths);
for t = 1:N
    xi=L*randn(n,nPaths);
    preturb=reshape(exp((mu-0.5*sigmas.^2)*dt.*ones(n,nPaths)+sigmas.*xi*sqrt(dt)),1,n,nPaths);
    S(t+1,:,:)=S(t,:,:).*preturb;
end
r=reshape((S(end,:,:)-S(1,:,:))./S(1,:,:),n,nPaths)';
end

