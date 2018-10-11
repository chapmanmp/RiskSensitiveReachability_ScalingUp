clear all
clc
close all
% script for value iteration with chain MDP

c = Chain_MDP_class(0.9,15);
tol = 1e-5;
maxIter = 1e3;

[V_Exp,Pol_Exp,err_Exp] = VI_Exp(c,tol,maxIter);
% 
% [V,Pol,err] = Value_Iteration(c,tol,maxIter);

%%% an index that controls the solver for inner optimization problem,
%%% either fmincon or CVX, CVX is slower in our example, but the
%%% conceptual, LP formulation can lead to something fast
index_opt = 'fmincon'; %'fmincon' or 'CVX'
Ny = 21;

% index_fast is an index that enables mex files for speed 1=ise mex file,
% 0=use matlab file
index_fast = 1; 
%% clear old mex files
if index_fast == 1 
% remove old files and build again
 rmdir('codegen','s');
 delete('P_interp_y_fast_mex.mexw64');

% Create a MEX configuration object
cfg = coder.config('mex');

% Create MEX file
arg_P = zeros(1,c.getNumStates);
arg_V = zeros(c.getNumStates,Ny);
arg_xi = zeros(c.getNumStates,1);
arg_yind = 1;
arg_y_set = zeros(1,Ny);
arg_state_ind = zeros(1,c.getNumStates);
arg_Ns_nz = 1;
codegen -config cfg P_interp_y_grad -args {arg_xi,arg_state_ind,arg_y_set,arg_V,arg_P,arg_yind,arg_Ns_nz}
disp('MEX generation complete!')

end
%% Do CVaR value iteration
Y_set_all = ones([c.Ns,1])*linspace(0,1,Ny);
 [V_CVaR,Pol_CVaR,err_CVaR] = VI_CVaR(c,Y_set_all,index_opt,index_fast,tol,maxIter);
 
%% plot stuff
V_plot= -V_CVaR;
Pol_plot = Pol_CVaR;
figure;
colormap = hsv(size(V_plot,2));

for i=1:10:size(V_plot,2)
% plot value function
figure(1)
plot(1:c.M, V_plot(1:c.M,i),'Color',colormap(i,:)); 
hold on
title('Value Function'); xlabel('state'); ylabel('value');

% plot policy
figure(2)
plot(1:c.M, Pol_plot(1:c.M,i),'Color',colormap(i,:)); 
title('Policy'); xlabel('state'); ylabel('action');
hold on
ylim([0,3]);

% plot the theoretical CVaR values
figure(3)
Y_span = linspace(0,1,Ny); % span of CVaR thresholds
r_theory = c.rewardCVaR(1:c.M,Y_span(i));

plot(1:c.M,r_theory,'Color',colormap(i,:));
hold on
ylim([0,4]);
title('Value Function'); xlabel('state'); ylabel('value');

% plot policy difference
pol_theory = find(r_theory == max(r_theory),1,'first');
figure(4);hold on;
pol_test = Pol_plot(1:c.M,i);
plot(1:c.M, (1:c.M < pol_theory)+1-pol_test','.','Color',colormap(i,:));
end

