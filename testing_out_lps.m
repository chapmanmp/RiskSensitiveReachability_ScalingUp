close all; clearvars; clc;

load('xis3over2_uis1_insidemaxExppond.mat');

%P = [0.0237; 0; 0; 0.5244; 0.3280; 0; 0; 0; 0; 0.1239]; % the P(i)'s that end up small are set to 0
%P = 1/10*ones(10,1);

%% solving lps one at a time

[tStar,bigexp] = argLargeLP_2( f_full, A, b, P, ls, nl, nd );
t1_smallLPs = tStar(1:nd); t2_smallLPs = tStar(nd+1:2*nd); 
t3_smallLPs = tStar(2*nd+1:3*nd); t4_smallLPs = tStar(3*nd+1:4*nd);

%% 1 larger LP - v2

bigA = []; for i = 1 : nl, bigA = blkdiag(bigA, A); end
bigb = repmat(b, nl, 1);
big_size = max(max(abs(b)))/2 + max(max(abs(A)))/2;
myfactor = max(P)/big_size;
cvx_solver mosek;
cvx_begin quiet
    cvx_precision best
    variables Z(nd,nl) t(nl*nd,1)   
    maximize( f_full' * t )  % (P/ls(i))' * t   
    subject to
    
        (myfactor*bigA) * vec(Z) + (myfactor*bigb) >= myfactor*vec( repmat(t', nl-1, 1) );
            
        P' * Z == ls';
         
        Z <= 1;
        Z >= 0;

cvx_end
t1_bigLP_compact = t(1:nd); t2_bigLP_compact = t(nd+1:2*nd); 
t3_bigLP_compact = t(2*nd+1:3*nd); t4_bigLP_compact = t(3*nd+1:4*nd);

%ti_bigLP_compact equals ti_smallLPs EXACTLY for each i! :)

%bigexp(1)-f_full(1:nd)'*t1_bigLP_compact ~ 10^(-23)
%bigexp(2)-f_full(nd+1:2*nd)'*t2_bigLP_compact = 0
%bigexp(3)-f_full(2*nd+1:3*nd)'*t3_bigLP_compact ~ 10^(-23)
%bigexp(4)-f_full(3*nd+1:4*nd)'*t4_bigLP_compact ~ 10^(-22)


%% 1 larger LP - more compact
almost_nl = 4;
bigA = []; for i = 1 : almost_nl, bigA = blkdiag(bigA, A); end
bigb = repmat(b, almost_nl, 1);
big_size = max(max(abs(b)))/2 + max(max(abs(A)))/2;
myfactor = max(P)/big_size;
cvx_solver mosek;
cvx_begin quiet
    cvx_precision best
    variables Z(nd,almost_nl) t(almost_nl*nd,1)   
    maximize( f_full(1:almost_nl*nd)' * t )  % (P/ls(i))' * t   
    subject to
    
        (myfactor*bigA) * vec(Z) + (myfactor*bigb) >= myfactor*vec( repmat(t', nl-1, 1) );
            
        P' * Z == ls(1:almost_nl)';
         
        Z <= 1;
        Z >= 0;

cvx_end
t1_bigLP_compact = t(1:nd); t2_bigLP_compact = t(nd+1:2*nd); t3_bigLP_compact = t(2*nd+1:3*nd);
%is exactly equal to loop over small LPs for almost_nl = 3

%for almost_nl = 4, t3_smallLPs-t3_bigLP_compact~10^-6 at entries 4 and 5, zero difference at other entries
    % ti_smallLPs-ti_bigLP_compact = 0 for i = 1, 2, 4
    % get the above result for cvx_precision default and cvx_precision best
    % P(4), P(5) ~ 0.4 whereas many other P(i) are 10^-4
    % Perhaps the inaccuracies are arising due to different weighting in
    % the objective function that is important when terms are summed -- but this doesn't arise in the smaller LP case
    % this might be fixed if we find P so that each >= 10^-2
    
%for almost_nl = 4, P >= 10^(-2), t3_smallLPs-t3_bigLP_compact~10^(-21), t4_smallLPs-t4_bigLP_compact~10^(-6)
%for almost_nl = 4, P([2,3,6:9]) = 0, t3_smallLPs-t3_bigLP_compact~10^(-6), t4_smallLPs-t4_bigLP_compact~10^(-7)
    % the errors are showing up where the probabilities are large
    
%for almost_nl = 4, P(i) = 1/10, t2_smallLPs-t2_bigLP_compact~10(-6) for entries 6-9
% t3_smallLPs-t3_bigLP_compact~10^(-6) for many entries -- t3* is on the order of 10^-6
% t4_smallLPs-t4_bigLP_compact~10^(-7) for many entries -- t4* is on the order of 10^-6
% because these errors are still happening for equivalent P(i) - I think
% that the problem is the different scales of (A,b) compared to (P,ls).
% Let's try to make the scales more similar. IT WORKED! Putting constraints
% on the same scales made ti_bigLP_compact equal to ti_smallLPs for each i,
% for original P :)

t4_bigLP_compact = t(3*nd+1:4*nd);

%% 1 larger LP

cvx_solver mosek;
cvx_begin quiet
    %cvx_precision best
    variables Z1(nd,1) Z3(nd,1) t1(nd,1) t3(nd,1)  
    maximize( f_full(1:nd)' * t1  + f_full(2*nd+1:3*nd)' * t3 )  % (P/ls(i))' * t   
    subject to           
            A*Z1 + b >= vec( repmat(t1', nl-1, 1) ); % [ti(1);...;ti(1);...;ti(nd);...;ti(nd)]
            %for i = 1 : nd,  As{i}*Z(i) + bs{i} >= t(i); end % one LMI per disturbance realization (eqv., per next state realization)
            A*Z3 + b >= vec( repmat(t3', nl-1, 1) );
            
            P' * Z1 == ls(1);
            P' * Z3 == ls(3);
            
            Z1 <= 1;
            Z1 >= 0;
            
            Z3 <= 1;
            Z3 >= 0; 
cvx_end
t1_bigLP_notcompact = t1; t3_bigLP_notcompact = t3;

