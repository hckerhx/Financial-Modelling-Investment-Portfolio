function [ x_optimal, otherInfo ] = diverse( mu, Q, targetRet, k )
% k: number of representatives
% otherInfo.z is n x n, shows which asset is represented by which

clear model model2

n = size(mu,1); 

% correlation coefficient 
rho = corrcov(Q);

% y1,y2...,z11,z12,z13,z14...z21,z22,z23,z24...
model.obj = [zeros(1,n), reshape(rho',1,[])];

% constriants
A=[ones(1,n), zeros(1, n^2);
    zeros(n,n), kron(eye(n),ones(1,n));
    -1*repmat(eye(n),n,1), eye(n^2)];
b=[k; ones(n,1);zeros(n^2,1)];

model.A = sparse(A);%repmat(kron(eye(r),[ones(1,1), zeros(1,r-1)]),r,1)]);
model.rhs = b;
model.sense = [repmat('=',1,n+1),repmat('<',1,n^2)];

% variable type
model.vtype = 'B';

model.modelsense = 'max';

params.outputflag = 0;

result = gurobi(model,params);

assert(size(result.status,2)==size('OPTIMAL',2),'Diverse Gurobi must converge')

x_result = result.x(1:n);

% upper bound and lower bound 
LB = -Inf(n,1);
LB(x_result<1)=0;

UB = Inf(n,1);
UB(x_result<1)=0;

% output 
x_optimal = MVO(mu, Q, targetRet, LB, UB);
otherInfo.z = reshape(result.x(n+1:n+n^2),n,n)';
otherInfo.totalCorr=result.objval;

end