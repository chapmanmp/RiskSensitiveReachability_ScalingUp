% Script for CVaR Value Iteration for small grid world domain

%% Grid-world parameters
eps = 0.05;         % noise parameter
target_row = 2;     % target location
target_col = 16;%60;
mdp = Gridworld_MDP_Pen_class('gridworld3.png',eps,target_row,target_col,'file',-2);

%% Do standard value iteration
maxIter = 1e3;
tol = 1e-5;
dis = 0.95;
[V_Exp,Pol_Exp,err_Exp] = VI_Exp(mdp,tol,maxIter,[],dis);
disp('Standard value iteration complete.');
% show values on image, and sample trajectory
start_row = 13;%50;
start_col = 16;%60;
im = mdp.val2image(-V_Exp(1:end-1));
figure;imagesc(im);colorbar;
[row,col] = mdp.getPath(Pol_Exp(1:end-1),start_row,start_col);
hold on;plot(col,row,'k-x','LineWidth',2,'MarkerSize',10);drawnow;

%% Do CVaR Value Iteration
% Choose solver: currently only CPlex is supported.
% 'linprog' - CPlex
index_opt = 'linprog'; 

% set VI parameters
Ny = 21;        % number of interpolation points for y
Y_set_all = ones([mdp.Ns,1])*[0, logspace(-2,0,Ny-1)]; 
maxIter = 40;
tol = 1e-3;
% Do CVaR VI
disp('Performing CVaR value iteration...');
[V_CVaR,Pol_CVaR,err_CVaR] = VI_CVaR(mdp,Y_set_all,index_opt,0,tol,maxIter,repmat(V_Exp,1,Ny),dis);
disp('done.');
%% Plot stuff
hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 1600 350])
kvals = [3,12,21];
for k = 1:numel(kvals)
    subplot(1,3,k); 
    im = mdp.val2image(-V_CVaR(1:end-1,kvals(k)));
    imagesc(im); colormap(parula(1e3));
    colorbar;
    [row,col] = mdp.getPath(Pol_CVaR(1:end-1,kvals(k)),start_row,start_col);
    hold on;plot(col,row,'k','LineWidth',2);
    title(['CVaR \alpha = ' num2str(round(Y_set_all(1,kvals(k)),2))],'fontsize',16);
end

