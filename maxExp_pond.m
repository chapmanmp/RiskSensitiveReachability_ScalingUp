%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: For each i, approximates max_R\inRiskEnvelope { E[ R*J_k+1(x_k+1,ls(i)*R) | x_k, ls(i), u_k ] }
%              Uses Chow 2015 linear interpolation method on confidence level
%              Uses change of variable, Z := y*R 
% INPUT: 
    % J_k+1 : optimal cost-to-go at time k+1, array
    % x : state at time k, real number
    % us : control options at time k, real number
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

function bigexp = maxExp_pond( J_kPLUS1, x, us, xs, ls, ws, P, dt, area_pond )

% # disturbance values  # confidence levels     # control options
nd = length(ws);        nl = length(ls);        nu = length(us);

f_full = zeros(nd,nl); for i = 1 : nl, f_full(:,i) = P/ls(i); end; f_full = vec(f_full); f_full = repmat(f_full, nu, 1);

nrows = nl*nd*(nl-1); bus = zeros(nrows,nu); fus = ones(nrows,nu); Aus = []; 
for j = 1 : nu
     % encodes linear interpolation of y*J_k+1( x_k+1, y ) versus y, given us(j) and x
    [ Au, bu ] = getLMI_pond( x, us(j), ws, xs, ls, J_kPLUS1, dt, area_pond );
    
    fu = max(P)/mean([max(max(abs(bu))), max(max(abs(Au)))]); % gets constraints on similar scales
    fu = 1;

    for i = 1 : nl, Aus = blkdiag(Aus, fu*Au); end
    
    bus(:,j) = repmat(fu*bu, nl, 1); fus(:,j) = fu*fus(:,j);

end
bus = vec(bus); fus = vec(fus);

%[tStar, bigexp] = argLargeLP_2( f_full, A, b, P, ls, nl, nd ); % column vector
%tStar = argLargeLP( f_full, A, b, P, ls, nl, nd ); % column vector
tStar = argLargeLP_us( f_full, Aus, bus, P, ls, nl, nd, nu, fus ); % column vector

%bigexp = zeros(nl,1);
bigexp = zeros(nl,nu);

for j = 1 : nu
    
    tStar_j = tStar( 1 + (j-1)*nl*nd : j*nl*nd ); % extract optimal arg for us(j)

    for i = 1 : nl,  bigexp(i,j) = (P/ls(i))' * tStar_j( (i-1)*nd + 1 : i*nd ); end % extract optimal value for us(j), ls(i)
    
end


