function [ x_optimal,Info ] = robust( mu, Q, targetRet, N, lambda, confidence )
% confidence = 0.9 means 90% confidence interval

clear model

[n,c] = size(mu);

% theta matrix
theta_matrix = diag(diag(Q))/N;

% epsilon
% here confidence = 1-alpha
epsilon = sqrt(chi2inv(confidence,n));

% set objective: x
% model.obj = [-mu' epsilon];
model.modelsense = 'min';

model.Q = sparse(lambda*blkdiag(Q,0));
model.obj=[-mu' epsilon];

% budget constraint 
A=[ones(1,n) 0];
b=1;
model.A   = sparse(A);
model.rhs = b;
model.sense = '=';

% second-order constraint
model.quadcon(1).Qc = sparse(blkdiag(theta_matrix,-1));
model.quadcon(1).q  = zeros(n+1,1);
model.quadcon(1).rhs = 0;
model.quadcon(1).sense = '=';

model.lb=[-1e4*ones(n,1);0];
params.OutputFlag = 0;

result = gurobi(model,params);
assert(size(result.status,2)==size('OPTIMAL',2),'Robust Gurobi must converge')
x_optimal = result.x(1:n);
Info.y=result.x(end);
Info.e=epsilon;
Info.theta=theta_matrix;
end