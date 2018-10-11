function [V,Pol,err] = VI_Exp(MDP,tol,maxIter,V0,dis)
% value iteration algorithm
% Inputs
% MDP : MDP class object
% tol : optional error tolerance (default 1e-3)
% maxIter : maximum iterations (default 1e3);
% V0 : initial value function (default 0)
% Outputs
% V : optimal value function
% Pol : optimal (greedy) policy

% Initialize
Ns = MDP.getNumStates;
if nargin < 5
    dis = 1;
end
if nargin < 4 || isempty(V0)
    V0 = zeros(Ns,1);
end
if nargin < 3
    maxIter = 1e3;
end
if nargin < 2
    tol = 1e-3;
end
V = V0;

Pol = zeros(Ns,1);

% Do Value Iteration
for i = 1:maxIter
    V_prev = V;
    for s = 1:Ns
        a = MDP.getActions(s);
        Q = -MDP.getReward(s,a) + dis*MDP.nextStateProb(s,a)*V_prev;
        [V(s),Pol(s)] = min(Q);
    end
    err = max(abs(V-V_prev))
    
    if err < tol
        disp(['Error threshold reached at iteration ' num2str(i)]);
        break;    
    end
end

end