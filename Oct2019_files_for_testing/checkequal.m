A = [pi 6; exp(1) 2]; b = [exp(2);8/pi]; P = [10^(-4); 1-10^(-4)]; ls = [.7 .4];

f_full = [P/ls(1); P/ls(2)];

cvx_solver mosek;
cvx_begin

    variables Z(2*2,1) t(2*2,1)
    
    maximize( f_full' * t )
    
    subject to
    
        A*Z(1:2) + b >= t(1:2);
        P'*Z(1:2) == ls(1);
        
        A*Z(3:4) + b >= t(3:4);
        P'*Z(3:4) == ls(2);
        
        Z <= 1;
        Z >= 0;
      
cvx_end
t1star = t(1:2); t2star = t(3:4);

val1_test = f_full(1:2)' * t1star; %equals val1 exactly for simple values of A, b, P
                                   % val1_test - val1 ~ 10^(-15) for less simple values

val2_test = f_full(3:4)' * t2star; %equals val2 exactly for simple values of A, b, P
                                   % val2_test - val2 ~ 10^(-15) for less simple values
%%

cvx_solver mosek;
cvx_begin

    variables Z1(2,1) t1(2,1)
    
    maximize( P' * t1 / ls(1) )
    
    subject to
    
        A*Z1 + b >= t1;
    
        Z1 <= 1;
        
        Z1 >= 0;
        
        P' * Z1 == ls(1);
        
cvx_end

t1star_test = t1;
val1 = cvx_optval;

%%

cvx_solver mosek;
cvx_begin

    variables Z2(2,1)  t2(2,1)
    
    maximize( P' * t2 / ls(2) )
    
    subject to
       
        A*Z2 + b >= t2;
 
        Z2 <= 1;
        
        Z2 >= 0;
      
        P' * Z2 == ls(2);

cvx_end
t2star_test = t2;
val2 = cvx_optval;

%%
K>> tStar(1:nd)

ans =

   1.0e-14 *

    0.2601
    0.2863
    0.3126
    0.3388
    0.3093
    0.3913
    0.4176
    0.4438
    0.4701
    0.4963

% summer code
maxExp_u1 =3.6509e-15


% scaling up code

BEFORE FIX
K>> maxExp_u1(1)
ans = 3.7590e-14

AFTER FIX
K>> maxExp_u1(1)
ans = 3.6509e-15

        
        