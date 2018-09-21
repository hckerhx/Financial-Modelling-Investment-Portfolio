function [ x_optimal, Info] = ...
    cvar( mu, Q, targetRet, nPaths, S0, T, confidence_level)

% scenarios: number of scenarios
% S0: 1 x n, initial price of n assets
% periods: periods of price paths to simulate, in this case 26 weeks (6
% months)
% 
% confidence=0.95 means 95% confidence interval; % 1-year window
N = 1; % Number of steps (half one week per time step)
n=size(S0,2);
% Find our correlation transformation matrix L
% L = chol(Q, 'lower');
% 
% % Number of simulated price paths
% % nPaths = 2000;
% 
% % Allocate Space for our Simulations:
% NumAssets = size(SO,2);
% 
% S_Assets = cell(1,NumAssets);
% 
% for i = 1:NumAssets
%     S_Assets{1,i} = [SO(i)*ones(1,nPaths); zeros(N, nPaths)];
% end
% 
% sigma = sqrt(diag(Q));
% 
% % Generate paths
% for k = 1:nPaths
%     for j = 1:N
%         for i = 1:NumAssets
%             % Convert the independent RVs to correlated RVs
%             coln = size(L, 2);
%             xi = L * randn(coln, 1); 
%             S_Assets{1,i}(j+1, k) = S_Assets{1,i}(j, k) * exp( ( mu(i) - 0.5 * sigma(i)^2 ) * dt ...
%                         + sigma(i) * sqrt(dt) * xi(i));
%         end
%     end
% end

[AssetRet S] = monteCarlo(S0,mu,Q,T,N,nPaths);
% gamma, z1,z2...zs, x1, x2...xs
f = [1, 1/((1-confidence_level)*nPaths)*ones(1,nPaths), zeros(1,n)];
A = [-ones(nPaths,1), -eye(nPaths), -AssetRet];
b = zeros(nPaths, 1);
Aeq = [0 zeros(1, nPaths) ones(1,n)];
beq = 1;
LB = [-Inf;zeros(nPaths+n,1)];
UB = Inf(nPaths+n+1,1);
[X,FVAL,EXITFLAG,OUTPUT] = linprog(f,A,b,Aeq,beq,LB,UB);
assert(EXITFLAG==1,'CVaR must converge')
x_optimal = X(end-n+1:end);
Info.r=AssetRet;
Info.S=S;
Info.cVaR=FVAL;
Info.VaR=X(1);

end



