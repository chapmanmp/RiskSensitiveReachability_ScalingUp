%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines the stage cost as an exponential of the signed distance w.r.t. constraint set
% INPUT: 
    % x: array of states (vector or matrix)
    % m: greater than or equal to 1, soft-max parameter
% OUTPUT: Stage cost
% AUTHOR: Margaret Chapman
% DATE: October 11, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = stage_cost_pond( x, m ) 

gx = signed_distance_pond( x );  % signed distance with respect to the constraint set

beta = 10^(-3);                  % scale down so LP CVX solver works for large m

c = beta * exp( m*gx );