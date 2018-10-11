%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Defines the signed distance function w.r.t. the constraint set, K = [0, 3)
% INPUT: State xk, or state trajectory (x0, x1, ..., xN)
% OUTPUT: Signed distance w.r.t. K = [0, 3)
% AUTHOR: Margaret Chapman
% DATE: September 5, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gx = signed_distance_pond( x )

gx = x - 3; % one-sided signed distance since we only care about overflowing and x will never be <0