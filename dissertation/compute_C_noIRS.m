function [Capacity_no_IRS] = compute_C_noIRS(Nr,H,P_max_linear,accuracy,sigma_lineal)

    
    %%% Singular value decomposition and acquisition of singular values

    [U4,D4,V4] = svd(H); 
    lam4 = diag(D4).^2; % Singular Values 


    %%% Water-filling power allocation

    wline = wfill(1./lam4, P_max_linear,accuracy);
    Pi = max(wline-1./lam4,0); % Tx power allocated to each Tx antenna
    power = zeros(1,size(V4,2));
    power(1:length(lam4)) = Pi; % optimal Tx power per Tx antenna + added zero to fit dimension of V;

    %%% Generation of the transmit covariance matrix

    Q = V4*diag(power)*V4'; % Tx covariance matrix

    %%% Channel capacity with IRS 

    Capacity_no_IRS = log2(det(eye(Nr)+1/sigma_lineal*H*Q*H'));

end



