% run Main_DynProgram to get Js
Js_new = Js;
clear Js
load('workspace_for_testing_Oct12.mat'); % 'ground truth' Js
Js_GROUNDTRUTH = Js;
clear Js

for i = 1 : length(Js_new)
    isequal(Js_new{i},Js_GROUNDTRUTH{i})
end

% 1 BIG LP (argLargeLP.m)
% Js_new{3}-Js_GROUNDTRUTH{3} ~ 10^(-14)
% Js_new{2}-Js_GROUNDTRUTH{2} ~ 10^(-7)
% Js_new{1}-Js_GROUNDTRUTH{1} ~ 10^(-7)

%Loop through many LPs (argLargeLP_2.m) P'*t/ls(i) or (P/ls(i))'*t
% Js_new{3}-Js_GROUNDTRUTH{3} ~ 10^(-14)
% Js_new{2}-Js_GROUNDTRUTH{2} ~ 10^(-10)
% Js_new{1}-Js_GROUNDTRUTH{1} ~ 10^(-10)


