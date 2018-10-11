function [S,A,Y,R,Rtotal] = Sample_CVaR_Traj_Perturbed(MDP,V,Y_set_all,alpha,T,S0,dis,MDP_P)
% Sample a T-step trajectory according to CVaR_alpha policy, from a perturbed MDP,
% with a given CVaR Value function
Ns = MDP.getNumStates;

% Set up optimization options and functions
options = cplexoptimset;
options.Display = 'off';

S = zeros(T,1);
S(1) = S0;
A = zeros(T,1);
R = zeros(T,1);
Y = zeros(T,1);
Y(1) = alpha;

for t = 1:T-1
    s = S(t);
    y_set = Y_set_all(s,:); % the discretization of y for each state
    a = MDP.getActions(s);
    objective_fn = zeros([length(a) 1]);
    xi_cell = cell(numel(a),1);
    for a_ind=1:numel(a)
        trans_prob = MDP.nextStateProb(s,a(a_ind));
        % intialize nz_prob_ind so it has a fixed length
        nz_prob_ind = zeros([1 length(trans_prob)]);
        % optimize only over non-zero probabilities
        ind_prob_pos = find(trans_prob > 0);
        nz_prob_ind(1:length(ind_prob_pos)) = ind_prob_pos; % rest is zero
        % number of non-zero elements
        Ns_nz = numel(ind_prob_pos);
        
        if Y(t) == 0
            % if y=0, we do not need to solve the inner optimization
            % problem because the corresponding result is V(x,0)
            objective_fn(a_ind) = max(V(ind_prob_pos));
            xi_cell{a_ind} = zeros(Ns,1);
        else
            lb = zeros([Ns_nz 1]);
            ub = @(y_in) ones([Ns_nz 1])/y_in;
            Aeq = trans_prob(ind_prob_pos);
            Beq = 1;
            state_ind = nz_prob_ind;
            state_ind(state_ind ==0) =[];
            Aeq_LP = [Aeq zeros([size(Aeq,1), Ns_nz])];
            Beq_LP = Beq;
            lb_LP = [lb;
                -1e6*ones([Ns_nz 1])];
            ub_LP = @(y_in) [ub(y_in);
                1e6*ones([Ns_nz 1])];
            A_LP = zeros([(length(y_set)-1)*Ns_nz 2*Ns_nz]);
            B_fn = cell(length(y_set)-1,1);
            %                 B_LP = zeros([(length(y_set)-1)*Ns_nz 1]);
            
            for n = 1:length(y_set)-1
                slope = (y_set(n+1)*V(state_ind,n+1)-y_set(n)*V(state_ind,n))/(y_set(n+1)-y_set(n));
                
                A_LP((n-1)*Ns_nz+1:n*Ns_nz,:) = [-diag(slope), eye(Ns_nz)];
                
                B_fn{n} = @(y_in) V(state_ind,n)*y_set(n)/y_in - slope*y_set(n)/y_in;
            end
            
            f_LP = [zeros([1 Ns_nz]) trans_prob(ind_prob_pos)];
            xi_0 = ones(Ns_nz, 1);
            B_LP = zeros([(length(y_set)-1)*Ns_nz 1]);
            for n = 1:length(y_set)-1
                B_LP((n-1)*Ns_nz+1:n*Ns_nz,1) = B_fn{n}(Y(t));
            end
            round_dec = 1000; % 3 decimal places for redundnat row elimination
            %get rid of dependent rows
            AB = round([A_LP,B_LP]*round_dec);
            [AB_unique,~,~] = unique(AB, 'rows', 'first');
            AB_unique = AB_unique/round_dec;
            A_LP_in = AB_unique(:,1:end-1);
            B_LP_in = AB_unique(:,end);
            tt =zeros([Ns_nz 1]);
            xi_t0 = [xi_0;tt];
            [xi_t, fval, exitflag, output] = cplexlp(-f_LP, A_LP_in, B_LP_in, Aeq_LP,Beq_LP, lb_LP,ub_LP(Y(t)), xi_t0,options );
            objective_fn(a_ind) = -fval;
            xi_cell{a_ind} = zeros(Ns,1);
            xi_cell{a_ind}(ind_prob_pos) = xi_t(1:Ns_nz);
        end
    end
    Q = -MDP.getReward(s,a) + dis*objective_fn(:);
    % calculate the greedy policy and value function update
    [~,A(t)] = min(Q);
    R(t) = MDP.getReward(S(t),A(t));
    xi_vec = xi_cell{A(t)};
    [row,col] = MDP.getRowCol(S(t));    % current row,col on original mdp
    S_P = MDP_P.getState(row,col);      % current state on perturbed mdp
    S_P_next = MDP_P.sampleNextState(S_P,A(t));
    [row,col] = MDP_P.getRowCol(S_P_next);      % next row,col on perturbed mdp
    S(t+1) = MDP.getState(row,col);             % next state on original mdp
    Y(t+1) = Y(t) * xi_vec(S(t+1));
end

Rtotal = sum((dis.^(1:T)).*R');

end


