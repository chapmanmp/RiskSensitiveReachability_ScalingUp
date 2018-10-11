function [V,Pol,err] = VI_CVaR(MDP,Y_set_all,index_opt,index_fast,tol,maxIter,V0,dis) 
% Interpolation based CVaR VI
% Instead of solving many small problem, solve a concatenated linear problem
% for all actions and confidence parameters at each iteration
Ns = MDP.getNumStates;
Ny = size(Y_set_all,2);
Na = MDP.getNumActions;
if nargin <8
    dis = 1;
end
if nargin <7
    V0 = zeros([Ns, Ny]);
end
if nargin < 6
    maxIter = 1e3;
end
if nargin < 5
    tol = 1e-3;
end
V = V0;
Pol = zeros(Ns,Ny);
warm_start = cell(Ns,Ny,Na);
f_cell = cell(Na,Ny-1);
A_cell = cell(Na,Ny-1);
B_cell = cell(Na,Ny-1);
Aeq_cell = cell(Na,Ny-1);
Beq_cell = cell(Na,Ny-1);
lb_cell = cell(Na,Ny-1);
ub_cell = cell(Na,Ny-1);
xi_t0_cell = cell(Na,Ny-1);

% Set up optimization options and functions
options = cplexoptimset;
options.Display = 'off';

% VI algorithm with Policy and Value function outputs
for i = 1: maxIter
    V_prev = V;
    for s = 1:Ns
       % I get rid of this line for better speed performance
       % disp(['iteration at state: ' num2str(Ns-s)]);
        y_set = Y_set_all(s,:); % the discretization of y for each state
        a = MDP.getActions(s);
        objective_fn = zeros([length(a) length(y_set)]);
        for a_ind=1:length(a)
            trans_prob = MDP.nextStateProb(s,a(a_ind));
            % intialize nz_prob_ind so it has a fixed length
            nz_prob_ind = zeros([1 length(trans_prob)]);
            % optimize only over non-zero probabilities
            ind_prob_pos = find(trans_prob > 0);
            nz_prob_ind(1:length(ind_prob_pos)) = ind_prob_pos; % rest is zero
            % number of non-zero elements
            Ns_nz = numel(ind_prob_pos);
            
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
            
            for n = 1:length(y_set)-1
                slope = (y_set(n+1)*V_prev(state_ind,n+1)-y_set(n)*V_prev(state_ind,n))/(y_set(n+1)-y_set(n));
                
                A_LP((n-1)*Ns_nz+1:n*Ns_nz,:) = [-diag(slope.*trans_prob(ind_prob_pos)'), eye(Ns_nz)];
                
                B_fn{n} = @(y_in) V_prev(state_ind,n)*y_set(n)/y_in - slope*y_set(n)/y_in;
                
            end
            
            f_LP = [zeros([1 Ns_nz]) ones([1 Ns_nz])];
            
            for y_ind = 1:length(y_set)
                
                if y_set(y_ind) == 0
                    % if y=0, we do not need to solve the inner optimization
                    % problem because the corresponding result is V(x,0)
                    objective_fn(a(a_ind),y_ind) = max(V_prev(ind_prob_pos,y_ind));
                    
                else
                    
                    if isempty(warm_start{s,y_ind,a_ind})
                        warm_start{s,y_ind,a_ind} = ones(Ns_nz, 1);
                    end
                    xi_0 = warm_start{s,y_ind,a_ind};
                    
                    B_LP = zeros([(length(y_set)-1)*Ns_nz 1]);
                    for n = 1:length(y_set)-1
                        B_LP((n-1)*Ns_nz+1:n*Ns_nz,1) = B_fn{n}(y_set(y_ind)).*trans_prob(ind_prob_pos)';
                    end
                    
                    t =zeros([Ns_nz 1]);
                    xi_t0 = [xi_0;t];
                    % Since cplexlp solves minimization problems and the problem
                    % is a maximization problem, negate the objective
                    f_cell{a_ind,y_ind-1} = -f_LP';
                    A_cell{a_ind,y_ind-1} = sparse(A_LP);
                    B_cell{a_ind,y_ind-1} = B_LP;
                    Aeq_cell{a_ind,y_ind-1} = sparse(Aeq_LP);
                    Beq_cell{a_ind,y_ind-1} = Beq_LP;
                    lb_cell{a_ind,y_ind-1} = lb_LP;
                    ub_cell{a_ind,y_ind-1} = ub_LP(y_set(y_ind));
                    xi_t0_cell{a_ind,y_ind-1} = xi_t0;
                end
            end
        end
        f_full = cell2mat(f_cell(:));
        B_full = cell2mat(B_cell(:));
        Beq_full = cell2mat(Beq_cell(:));
        lb_full = cell2mat(lb_cell(:));
        ub_full = cell2mat(ub_cell(:));
        xi_t0_full = cell2mat(xi_t0_cell(:));
        A_concat = A_cell(:);
        A_full = sparse(0,0);
        for nn = 1:length(A_concat)
            A_full = blkdiag(A_full,A_concat{nn});
        end
        Aeq_concat = Aeq_cell(:);
        Aeq_full = sparse(0,0);
        for nn = 1:length(Aeq_concat)
            Aeq_full = blkdiag(Aeq_full,Aeq_concat{nn});
        end
        % solve optimization problem
        [xi_t, fval, exitflag, output] = cplexlp(f_full, A_full, B_full, Aeq_full,Beq_full, lb_full, ub_full, xi_t0_full, options );
        % distribute result to states
        cellsz = cellfun(@numel,xi_t0_cell,'uni',false);
        xi_cell = reshape(mat2cell(xi_t,cell2mat(cellsz(:))),size(xi_t0_cell));
        objective_fn(:,2:end) = cellfun(@(x,y) -x'*y, f_cell, xi_cell);
        for y_ind = 1:length(y_set)
            % form Q function
            Q = -MDP.getReward(s,a) + dis*objective_fn(:,y_ind);
            
            % calculate the greedy policy and value function update
            [V(s,y_ind),Pol(s,y_ind)] = min(Q);
            
            ind_min = find(Q==V(s,y_ind));
            
            if length(ind_min)>1
                % disp('more than one optimal policy')
                Pol(s,y_ind) = max(ind_min);
            else
                Pol(s,y_ind) = ind_min;
            end
        end
    end
    err = max(reshape(abs(V-V_prev),[size(V,1)*size(V,2),1]))
    if err < tol
        disp(['Error threshold reached at iteration ' num2str(i)]);
        break;
    end
    
    
end


