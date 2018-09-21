function  x_optimal = MVO(mu, Q, targetRet,lb,ub)

% Use this function to construct your MVO portfolio subject to the
% target return, with short-selling allowed.
%
% You may use quadprog, Gurobi, or any other optimizer you are familiar
% with. Just be sure to comment on your code to (briefly) explain your
% procedure.

% Find the total number of assets
n = size(Q,1);
if nargin<=3
    lb=-Inf(n,1);
    ub=Inf(n,1);
elseif  nargin<=4
    ub=Inf(n,1);
end

% *************** WRITE YOUR CODE HERE ***************
%----------------------------------------------------------------------
% Quadprog
% must be greater than target return
A = -mu'; % negative because in std form Ax<b
b = -targetRet;
% budget constraint
Aeq = ones(1,n);
beq = 1;

[X FVAL EXITFLAG]=quadprog(Q,[],A,b,Aeq,beq,lb,ub);
assert(EXITFLAG==1,'MVO must converge')

x_optimal = X;


% Gurobi should return the same result
% must be greater than target return
% A = -mu'; % negative because in std form Ax<b
% b = -targetRet;
% % budget constraint
% Aeq = ones(1,n);
% beq = 1;
%
% clear model;
% intcon='C';
% rowcon=['<','='];
%
% model.Q = sparse(Q);
% model.A =sparse([A;Aeq]);
% model.obj=zeros(1,n);
% model.rhs = [b;beq];
% model.sense = char(rowcon);
% model.lb=-Inf*ones(1,n);
% % no need to set upper bound, default infinite
%
% gurobi_write(model, 'mvo.lp');
%
% model.vtype = char(intcon);
%
%
% results  = gurobi(model);
%
%
% x_quad = results.x(1:n);

% assert(sum(x_quad-x_optimal)<0.001,'Two methods should return the same weights')
%----------------------------------------------------------------------

end