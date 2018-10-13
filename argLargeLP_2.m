%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Returns optimal argument to compute for each i,
%                   max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k,ls(i), u_k ]
%              Uses Chow 2015 linear interpolation method on confidence level
%              Uses change of variable, Z := y*R
% INPUT:
    % f_full = [P/ls(1); P/ls(2); ...; P/ls(end)] column vec
    % A (matrix), b (col vector) : encode linear interpolation of y*J_k+1( x_k+1, y ) vs. y given x and u
    % P(i): probability that w_k = ws(i)
    % ls(i): ith confidence level
    % nl = number of confidence levels
    % nd = number of disturbance values, length of P
% OUTPUT:
    % [t1; t2; ...] (col vector) that maximizes 1/ls(1)P'*t1 + 1/ls(2)P'*t2 + ...
        % subject to constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tStar,bigexp] = argLargeLP_2( f_full, A, b, P, ls, nl, nd )
tStar = []; bigexp = zeros(nl,1);
for i = 1 : nl
start_i = (i-1)*nd + 1; end_i = i*nd;
cvx_solver mosek;
cvx_begin quiet
    %cvx_precision best
    variables Z(nd,1) t(nd,1)   
    maximize( (f_full(start_i:end_i))' * t )  % (P/ls(i))' * t   
    subject to           
            A*Z + b >= vec( repmat(t', nl-1, 1) ); % [ti(1);...;ti(1);...;ti(nd);...;ti(nd)]
            %for i = 1 : nd,  As{i}*Z(i) + bs{i} >= t(i); end % one LMI per disturbance realization (eqv., per next state realization)
            P' * Z == ls(i);
            
            Z <= 1;
            
            Z >= 0;    
cvx_end  
tStar = [tStar; t];
bigexp(i) = cvx_optval;

if isinf(cvx_optval) || isnan(cvx_optval) || strcmpi(cvx_status, 'Inaccurate/Solved')
    display(cvx_optval);
    error('maxExp.m: cvx not solved.'); 
end

end
end