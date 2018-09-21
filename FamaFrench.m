function result = FamaFrench( periodReturns,periodFactRet )
n=size(periodReturns,2);
params = zeros(4, n);
numweeks = size(periodReturns,1);
%     each column is residuals for an asset
epsilons = zeros(numweeks,n);
for i=1:n
    % For each asset, we have 4 factor loadings:
    % a(active return)
    % beta Market
    % beta SMB
    % beta HML
    [params(:,i),interval,epsilons(:,i)]=...
        regress(periodReturns(:,i),[ones(numweeks,1),periodFactRet]);
end

A = params(1,:)';          % n x 1 vector of alphas
V = params(2:end,:);       % m x n matrix of betas

%     geo-mean returns a row vector of each column's mean, need transpose

f_bar = (Geomean(periodFactRet+1)-1)';  % m x 1 vector of factor expected returns
F = cov(periodFactRet);          % m x m factor covariance matrix

epsilon = epsilons;    % Regression residuals
D = cov(epsilon);          % Diagonal n x n matrix of residual variance
%     Note: D is not diagonal but OK, (Fama-french model does not respect
%     ideal environment

mu = A+V'*f_bar;         % n x 1 vector of asset exp. returns
Q_1  = V'*F*V + D;         % n x n asset covariance matrix
%     Note: Q_1 is not completely symmetric, truncation error
Q = cov(periodReturns);
assert(max(max(Q-Q_1))<0.001,'Q should be same with or without using factor model')
result.mu=mu;
result.Q=Q;
end

