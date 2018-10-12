% run Main_DynProgram to get Js
Js_new = Js;
clear Js
load('workspace_for_testing_Oct11.mat'); % 'ground truth' Js
Js_old = Js;
clear Js

for i = 1 : length(Js_new)
    isequal(Js_new{i},Js_old{i})
end