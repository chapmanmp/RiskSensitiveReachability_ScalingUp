function [val, grad] = P_interp_y_grad(xi,state_ind,y_set,V,P,y_ind,Ns_nz) %#codegen
% from state (s,y) to next state x interpolation function, also
% divided by y. grad is the gradient w.r.t. xi.
grad = zeros([Ns_nz 1]);
P_short = P(state_ind(1:Ns_nz));
% Note: y_ind here makes sure that y is non-zero !

output = zeros([Ns_nz 1]);
% output is an interpolation of yV(x,y), i.e. I[yV], and divides by y
xi_short = xi(1:Ns_nz);
for x = 1:numel(xi_short) % loop over states with non-zero prob, length = Ns_nz
    if xi_short(x) == 0
        output(x) = 0;
        grad(x) = 0;
    else  % xi(x) \~= 0 and y~=0
        ind = find(y_set(y_ind)*xi_short(x)-y_set>=0,1,'last');
        % adding a tol parameter is to account for the
        % numerical error in cut-off
        if isempty(ind)
            ind = 1; % This is due to the numerical error in thresholding
        elseif ind(1) == length(y_set)
        % to make sure ind will not crash if y_set(y_ind)*xi_short(x)=1
            ind = ind -1;
        end
        
        slope = (y_set(ind+1)*V(state_ind(x),ind+1)-y_set(ind)*V(state_ind(x),ind))/(y_set(ind+1)-y_set(ind));
        output(x) = (y_set(ind)*V(state_ind(x),ind) + slope*(y_set(y_ind)*xi_short(x)-y_set(ind)))/y_set(y_ind); % interpolation function  divide by y
        grad(x) = slope;
    end
end
val = -P_short*output; % minus accounts for maximization
grad = -P_short'.*grad;
end