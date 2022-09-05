
            %%%% Sum-rate optimization for IRS-aided MIMO %%%%

    clc;
    clear;

 %%% Definition of all the variables
 
    M = 40;
    f = 10e6; 
    c = 3e8;
    lambda = c/f;
    theta_AoA = 0;
    theta_AoD = 0;
    Nt = 4; %% Number of Tx-IRS streams
    Nr = 4; %% Number of IRS-Rx streams
    a = 0; %% Lower-bound reflection coefficient phase in degrees
    b = 360; %% Upper-bound reflection coefficient phase in degrees
    H_hat = 10;
    dp_hat = 2;
    dh_hat = 2;
    do = 1;
    alpha_IR = 2.8;
    alpha_TI = 2.2;
    alpha_d = 3.5;
    dD_hat = 600;
    Kd = 0;
    beta_o = -30;
    beta_o_lineal = 10^(-30/10);
    dA = lambda/2;
    P_max_dBm = 30; % in dBm
    P_max_dB = P_max_dBm-30; % in dB
    P_max_linear = 10^(P_max_dB/10);
    accuracy = 1E-3; % accuracy of water-filing algorithm
    sigma_dBm = -90;
    sigma_dB = sigma_dBm-30;
    sigma_lineal = 10^(sigma_dB/10);
 
%%% Definition of all the channels and inizialization

    
[R] = IRS_Receiver_channel(H_hat,Nr,dp_hat,dh_hat,do,alpha_IR,M,beta_o_lineal);
[T] = Transmitter_IRS_channel(Nt,M,dp_hat,dh_hat,dD_hat,do,alpha_TI,beta_o_lineal);
[H] = direct_channel(H_hat,dD_hat,Nr,Nt,lambda,theta_AoA,theta_AoD,Kd,do,beta_o_lineal,dA,alpha_d);


[optR,optQ,totalQ,totalC,totalR,Q,H_tilde,reflection_coefficients,optCapacity,Capacity_IRS] = computeQ_R(Nr,Nt,M,R,T,H,P_max_linear,accuracy,sigma_lineal,a,b);

Capacity_IRS_random_phase = log2(det(eye(Nr)+1/sigma_lineal*H_tilde*Q*H_tilde'));  %%% Capacity for a given set of reflection coefficients with random phase

[Capacity_no_IRS] = compute_C_noIRS(Nr,H,P_max_linear,accuracy,sigma_lineal);  %%% This function gives the channel capacity of the system without IRS

%%% Up to this point: Optimal transmit covariance matrix and reflection
%%% coefficients initialized (OptR and OptQ)


cond = true;
NumIterations = 1;
CapacityVec = [];
IterationsVec = [];

capacity_prev = 0;

optimal_capacity_11 = [];

while(cond)  %%% This process will be executed until two adjacent capacity values are the same
   
%%% Computing Am and Bm - Subproblem 2

[Am,Bm,Final_matrix,Am1,Bm1,NZ_eigenvalues,optimal_solution_1,optimal_capacity_1] =  compute_Am_Bm(M,R,T,H,optR,optQ,Nr,Nt,sigma_lineal);

%%% Computing optimized Q - Subproblem 1

[Q1,reflection_coefficients_1,H_tilde_1] = compute_final_Q(Nr,Nt,M,R,T,H,P_max_linear,accuracy,optimal_solution_1,sigma_lineal);

Capacity_IRS_1 = log2(det(eye(Nr)+1/sigma_lineal*H_tilde_1*Q1*H_tilde_1'));


optQ = Q1;
optR = reflection_coefficients_1;


    if abs(Capacity_IRS_1 - capacity_prev) < 10e-4  %%% Iteratively comparing two adjacent capacity values
            cond = false;  
    end
    
capacity_prev = Capacity_IRS_1;


NumIterations = NumIterations + 1;


end
 
%%% disp(NumIterations);
%%% disp(Capacity_IRS_1);
%%% disp(Capacity_IRS);
%%% disp(abs(Capacity_IRS_1 - capacity_prev));
%%% disp(Capacity_no_IRS);
%%% disp(Capacity_IRS_random_phase);


