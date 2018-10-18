%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Set-up script for pond example, time-invariant finite probability distribution
    % wk : average surface runoff rate on [k,k+1), [ft^3/s]
% AUTHOR: Margaret Chapman
% DATE: September 5, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_lb = 0; K_ub = 3;                     % Constraint set bounds [ft], K = (0, 3ft)

dx = 1/2;                               % State discretization [ft]

xs = K_lb : dx : K_ub + 1;              % Discretized states [ft]
nx = length(xs);
         
ls = [0.95; 0.7; 0.3; 0.05];            % Discretized confidence levels

[ X, ~ ] = meshgrid( xs, ls' );

dt = 300;                               % Duration of [k, k+1) [sec], 5min = 300sec

T = 4*3600;                             % Design storm length [sec], 4h = 4h*3600sec/h

%N = T/dt;                              % Time horizon: {0, 1, 2, ..., N} = {0, 5min, 10min, ..., 240min} = {0, 300sec, 600sec, ..., 14400sec}
N = 3;

surface_runoff_stats;                   % Provides possible values of surface runoff (ws) & moments
                                        % ws(i): ith possible value of wk (ft^3/s)

P = getProbDist(ws, Mymean, Myvariance, Myskewness); % P(i): probability that wk = ws(i)

m = 10;                                 % soft-max parameter

area_pond = 28292;                      % approx. surface area, pond 1 (south) [ft^2]

