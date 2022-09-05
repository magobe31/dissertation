
%%%% Resolution of Subproblem 1: Optimization of Q with given 
%%%% reflection coefficients  

function [Q1,reflection_coefficients_1,H_tilde_1] = compute_final_Q(Nr,Nt,M,R,T,H,P_max_linear,accuracy,optimal_solution_1,sigma_lineal)

  reflection_coefficients_1 = optimal_solution_1;
  
  E = zeros(Nr,Nt);

    for n = 1:M
        E1 = reflection_coefficients_1(n)*R(:,n)*T(n,:);
        E = E + E1; 
    end
    

    %%% Generation of the channel with IRS

    H_tilde_1 = H + E;

    %%% Singular value decomposition and acquisition of singular values

    [U3,D3,V3] = svd(H_tilde_1); 
    lam1 = diag(D3).^2; %Singular Values 

     %%% Water-filling power allocation

    wline = wfill(1./lam1, P_max_linear,accuracy);
    Pi1 = max(wline-1./lam1,0); % Tx power allocated to each Tx antenna
    power1 = zeros(1,size(V3,2));
    power1(1:length(lam1)) = Pi1; % optimal Tx power per Tx antenna + added zero to fit dimension of V;
    
     %%% Generation of the transmit covariance matrix

    Q1 = V3*diag(power1)*V3'; % Tx covariance matrix
     
end