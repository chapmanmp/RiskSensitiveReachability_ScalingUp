%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Performs the CVaR Bellman backwards recursion, uk \in {0,1}
% INPUT: 
    % J_k+1 : optimal cost-to-go at time k+1, array
    % X : row of states repeated [ xs; xs; ... ], array
    % L : column of confidence levels repeated [ ls ls ... ], array
    % ws(i): ith possible value of wk
    % P(i): probability that wk = ws(i)
    % m : soft-max parameter for stage_cost_pond.m
    % dt : duration of [k,k+1) interval
    % area_pond : approx. surface area of pond
% OUTPUT: 
    % J_k : optimal cost-to-go starting at time k, array
    % mu_k : optimal controller at time k, array
% AUTHOR: Margaret Chapman
% DATE: September 5, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ J_k, mu_k ] = CVaR_Bellman_Backup_pond( J_kPLUS1, xs, nx, ls, ws, P, m, dt, area_pond )

J_k = J_kPLUS1; mu_k = J_kPLUS1; % initialization

for i = 1 : nx      % <--x's change along columns of J_k, X, L-->
    
    x = xs(i);
    
    us = [0; 1];    % possible control actions at state x
    
    maxExp_us = maxExp_pond( J_kPLUS1, x, us, xs, ls, ws, P, dt, area_pond ); 
    %maxExp_us(i,j) = maxExp, given state x, confidence level ls(i), and control us(j)
    
    %maxExp_u1 = maxExp_pond( J_kPLUS1, x, us(1), xs, ls, ws, P, dt, area_pond ); 
    %maxExp_u2 = maxExp_pond( J_kPLUS1, x, us(2), xs, ls, ws, P, dt, area_pond );
    %[ optExp, optInd ] = min( [maxExp_u1, maxExp_u2], [], 2 ); % col vector, 1 entry per confidence level
   
    [ optExp, optInd ] = min( maxExp_us, [], 2 ); % col vector, 1 entry per confidence level
    
    J_k(:,i) = stage_cost_pond(x,m) + optExp;
    
    mu_k(:,i) = us(optInd); % need to check
    
end