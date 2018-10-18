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
    % nu = number of control options 
% OUTPUT:
    % [t1; t2; ...] (col vector) that maximizes 1/ls(1)P'*t1 + 1/ls(2)P'*t2 + ...
        % subject to constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tStar = argLargeLP_us( f_full, bigA, bigb, P, ls, nl, nd, nu, fus )

%big_size = max(max(abs(bigb)))/2 + max(max(abs(bigA)))/2; myfactor = max(P)/big_size; % gets constraints on similar scales

cvx_solver mosek;
cvx_begin quiet
    cvx_precision best
    
    variables Z(nd,nl*nu) t(nl*nd*nu,1)
    
    maximize( f_full' * t )  % (P/ls(i))' * t   
    subject to
        
        %(myfactor*bigA) * vec(Z) + (myfactor*bigb) >= myfactor*vec( repmat(t', nl-1, 1) );
        
        bigA * vec(Z) + bigb >= fus.*vec( repmat(t', nl-1, 1) ); %bigA, bigb have fus incroporated
        P' * Z == repmat( ls', 1, nu ); %ls is a column vector
        Z <= 1;
        Z >= 0;

cvx_end

if isinf(cvx_optval) || isnan(cvx_optval) || ~strcmpi(cvx_status, 'Solved')
    error('maxExp.m: cvx not solved.'); 
end

tStar = t;


