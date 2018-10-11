function mdp_per = perturb_mdp(mdp, ep, target_row, target_col)
mdp_p = copy(mdp);
ob_ind = find(mdp_p.img == 0);
for i = 1:numel(ob_ind)
    if rand < ep
        [row,col] = ind2sub(size(mdp_p.img),ob_ind(i));
        r = rand;
        if r < 0.25
            ind = sub2ind(size(mdp_p.img),max(row-1,1),col);
            mdp_p.img(ob_ind(i)) = 255;
            mdp_p.img(ind) = 0;
        elseif r < 0.5
            ind = sub2ind(size(mdp_p.img),min(row+1,size(mdp_p.img,1)),col);
            mdp_p.img(ob_ind(i)) = 255;
            mdp_p.img(ind) = 0;
        elseif r < 0.75
            ind = sub2ind(size(mdp_p.img),row,max(col-1,1));
            mdp_p.img(ob_ind(i)) = 255;
            mdp_p.img(ind) = 0;
        else
            ind = sub2ind(size(mdp_p.img),row,min(col+1,size(mdp_p.img,2)));
            mdp_p.img(ob_ind(i)) = 255;
            mdp_p.img(ind) = 0;
        end
    end
end
mdp_per = Gridworld_MDP_Pen_class( mdp_p.img,mdp.eps,target_row, target_col,'matlab');
end