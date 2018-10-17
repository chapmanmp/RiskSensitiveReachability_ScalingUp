% run Main_DynProgram to get Js
Js_new = Js;
clear Js
load('workspace_for_testing_Oct15.mat'); % 'ground truth' Js
Js_GROUNDTRUTH = Js;
clear Js

for i = 1 : length(Js_new)
    isequal(Js_new{i},Js_GROUNDTRUTH{i})
end

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



