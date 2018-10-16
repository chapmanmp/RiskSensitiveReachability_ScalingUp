%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: For each i, approximates max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), u_k ] }
%              Uses Chow 2015 linear interpolation method on confidence level
%              Uses change of variable, Z := y*R 
% INPUT: 
    % J_k+1 : optimal cost-to-go at time k+1, array
    % x : state at time k, real number
    % u : control at time k, real number
    % xs : x values, row vector
    % ls : confidence levels, row vector
    % ws(i): ith possible value of w_k
    % P(i): probability that w_k = ws(i)
    % dt : duration of [k,k+1) interval
    % area_pond : approx. surface area of pond
% OUTPUT: bigexp(i) ~= max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), u_k ] }
% AUTHOR: Margaret Chapman
% DATE: October 12, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bigexp = maxExp_pond( J_kPLUS1, x, u, xs, ls, ws, P, dt, area_pond )

% # disturbance values  # confidence levels
nd = length(ws);        nl = length(ls);

f_full = []; for i = 1 : nl, f_full = [f_full; P/ls(i)]; end

[ A, b ] = getLMI_pond( x, u, ws, xs, ls, J_kPLUS1, dt, area_pond );
% encodes linear interpolation of y*J_k+1( x_k+1, y ) versus y, given u and x

%[tStar, bigexp] = argLargeLP_2( f_full, A, b, P, ls, nl, nd ); % column vector
tStar = argLargeLP( f_full, A, b, P, ls, nl, nd ); % column vector

bigexp = zeros(nl,1);

for i = 1 : nl 
 
    start_i = (i-1)*nd + 1; end_i = i*nd;

    bigexp(i) = f_full(start_i:end_i)' * tStar(start_i:end_i);

end
