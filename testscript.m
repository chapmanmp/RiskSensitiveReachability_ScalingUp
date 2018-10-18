% run Main_DynProgram to get Js
Js_new = Js;
clear Js
load('results_acc2019\dyn_prog_results_sept112018.mat'); % 'ground truth' Js
Js_GROUNDTRUTH = Js;
clear Js

max(max(abs(Js_new{1}-Js_GROUNDTRUTH{1})))

% 1 BIG LP (argLargeLP_us.m) for all ls and us, no scaling
% on ACC example (66 states, 9 confidence levels, 48 time points, etc.)
% max(max(abs(Js_new{1}-Js_GROUNDTRUTH{1}))) = 0.1929



%---------------------------------------------
% 1 BIG LP (argLargeLP_us.m) for all ls and us, less bulk scaling
% Js_new{3}-Js_GROUNDTRUTH{3} ~ 1.0e-14
% Js_new{2}-Js_GROUNDTRUTH{2} ~ 1.0e-13
% Js_new{1}-Js_GROUNDTRUTH{1} ~ 1.0e-09 - showing up at xs(5), ls(2) - I'm not sure why this is happening

% 1 BIG LP (argLargeLP_us.m) for all ls and us, bulk scaling
% Js_new{3}-Js_GROUNDTRUTH{3} ~ 1.0e-14
% Js_new{2}-Js_GROUNDTRUTH{2} ~ 1.0e-12
% Js_new{1}-Js_GROUNDTRUTH{1} ~ 1.0e-09

% above results indicate that bulk scaling is okay
%---------------------------------------------

% 1 BIG LP (argLargeLP.m) for all ls
% Js_new{3}-Js_GROUNDTRUTH{3} ~ 1.0e-15
% Js_new{2}-Js_GROUNDTRUTH{2} ~ 1.0e-13
% Js_new{1}-Js_GROUNDTRUTH{1} ~ 1.0e-12

% Loop through many LPs (argLargeLP_2.m) (P/ls(i))'*t with bigexp(i) taken directly from cvx_optval
% Js_new{i} = Js_GROUNDTRUTH{i} exactly



