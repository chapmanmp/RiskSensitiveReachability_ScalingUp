% run a perturbation experiment script
k_vec = [Y_set_all(1,3),Y_set_all(1,7),Y_set_all(1,12),Y_set_all(1,21)];
t_row = 2;
t_col = 16;
num_per = 20;
num_runs = 20;
res = zeros(num_per, num_runs, numel(k_vec));
figure;imagesc(mdp.img);
for i = 1:num_per
    disp(['perturbation ' num2str(i)]);
    mdp_p = perturb_mdp(mdp, 0.5, target_row, target_col);
    for j = 1:num_runs 
        for k = 1:numel(k_vec)
            [S,A,Y,R,Rtotal] = Sample_CVaR_Traj_Perturbed(mdp,V_CVaR,Y_set_all,k_vec(k),100,mdp.getState(start_row,start_col),dis,mdp_p);
            res(i,j,k) = Rtotal;
            [rowC,colC] = mdp.getRowCol(S);
            if k == 1 
                hold on;plot(colC,rowC,'r');
            elseif k == 2
                hold on;plot(colC,rowC,'k');
            elseif k == 3
                hold on;plot(colC,rowC,'c');                
            else
                hold on;plot(colC,rowC,'m');
            end
            drawnow;
        end
    end
end
%%
r_dist = zeros(numel(k_vec),num_per*num_runs);
for k = 1:numel(k_vec)
    r_vec = squeeze(res(:,:,k));
    r_dist(k,:) = r_vec(:);
end