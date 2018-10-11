classdef Chain_MDP_class < Finite_MDP_class
    % Chain MDP chain MDP with state space $X = \{ 1,2,3,\dots,M \}$. 
    % For each state $x\in X$ there are two actions $u_1$ and $u_2$. 
    % Under $u_1$ the agent receives reward $r(x) = x$, and the 
    % process terminates (state 0). Under $u_2$ the state transitions to 
    % $x+1$ w.p. $p$ and zero reward, and w.p. $1-p$ the process terminates, 
    % also with zero reward.
    properties
        p = [];     % transition probability
        M = 1;
    end
    methods 
        function obj = Chain_MDP_class(p,M)
            Ns = M + 1;
            Na = 2;
            % build reward vector
            R = zeros(Ns,Na);
            R(1:M,1) = (1:M);
            % build transition matrix
            P1 = zeros(Ns,Ns);
            P1(:,M+1) = 1;
            P2 = zeros(Ns,Ns);
            for i = 1:M-1
                P2(i,i+1) = p;
                P2(i,M+1) = 1-p;
            end
            P2(M,M+1) = 1;
            P2(M+1,M+1) = 1;
            P = cat(3,P1,P2);
            obj@Finite_MDP_class(P,R);
            obj.p = p;
            obj.M = M;
        end
        function [r,prob] = rewardDist(obj,pol_th)
            % calculate the reward distribution for a threshold policy that
            % selects action 2 for s < pol_th and action 1 when s >= pol_th
            r = 0:obj.M;
            prob = zeros(1,obj.M+1);
            prob(pol_th+1) = obj.p.^(pol_th-1); % success probability
            prob(1) = obj.p.^(pol_th-1); % else zero reward
        end
        function [rmean] = rewardMean(obj,pol_th)
            % calculate the reward mean for a threshold policy that
            % selects action 2 for s < pol_th and action 1 when s >= pol_th
            % if pol_th is a vector then a vector of means is returned
            rmean = zeros(size(pol_th));
            for i = 1:length(pol_th)
%                 [r,prob] = obj.rewardDist(pol_th(i));
%                 rmean(i) = prob*r';
                rmean(i) = obj.p.^(pol_th(i)-1) * pol_th(i);
            end
        end
        function [rcvar] = rewardCVaR(obj,pol_th,alpha)
            % calculate the reward CVaR_alpha for a threshold policy that
            % selects action 2 for s < pol_th and action 1 when s >= pol_th
            % if pol_th is a vector then a vector of CVaRs is returned
            rcvar = zeros(size(pol_th));
            for i = 1:length(pol_th)
                ps = obj.p.^(pol_th(i)-1);
                rcvar(i) = ((alpha-1+ps)/alpha) * pol_th(i);
            end
        end
        function [val] = valueMean(obj)
            % calculate optimal value function
            rmean = obj.rewardMean(1:obj.M);
            pol = find(rmean == max(rmean),1,'first');
            val = 1:obj.M;
            for i = 1:pol
                val(i) = obj.p.^(pol-i) * pol;
            end
        end
        function [val] = valueCVaR(obj,alpha)
            % calculate optimal CVaR_alpha value function
            rmean = obj.rewardCVaR(1:obj.M,alpha);
            pol = find(rmean == max(rmean),1,'first');
            val = 1:obj.M;
            for i = 1:pol
                ps = obj.p.^(pol-i);
                val(i) = ((alpha-1+ps)/alpha) * pol;
            end
        end
    end
end