%Scaling experiments

%# time points, N = 48
%# confidence levels, nl = 9
%# controls, nu = 2
%# disturbances, nd = 10

%solving N*nx LPs, where each LP has 2*nd*nl*nu variables

%# states,   nx = 7   (dx = 1  ): time = 1  min 17 sec
 %           nx = 14  (dx = 1/2): time = 2  min 39 sec
  %          nx = 27  (dx = 1/4): time = 5  min 1 sec
   %         nx = 53  (dx = 1/8): time = 9  min 33 sec
    %        nx = 66 (dx = 1/10): time = 12 min 4.70 sec
     %      nx = 105 (dx = 1/16): time = 23 min
      %     nx = 209 (dx = 1/32): time = 39 min 3.70 sec
       %    nx = 417 (dx = 1/64): time = 1 hr 20 min
        %  nx = 833 (dx = 1/128): time = 2 hr 39 min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_states = [7, 14, 27, 53, 66, 105, 209, 417, 833];

compute_time = [1 + 17/60, 2 + 39/60, 5 + 1/60, 9 + 33/60, 12 + 4.70/60, 23, 39 + 3.70/60, 80, 120 + 39]; % min

figure; FigureSettings; plot(num_states, compute_time, ':^k');

ylabel('Computation time (min)'); legend('uniform grid'); xlabel('Number of states'); grid on;