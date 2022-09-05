
%%%% Initialization process: Optimizing Q given a set of reflection
%%%% coefficients with random phase and selecting Q and reflection coeff.
%%%% giving the maximum achievable capacity as starting point

function [optR,optQ,totalQ,totalC,totalR,Q,H_tilde,reflection_coefficients,optCapacity,Capacity_IRS] = computeQ_R(Nr,Nt,M,R,T,H,P_max_linear,accuracy,sigma_lineal,a,b)

L = 100;
totalQ = [];
totalC = [];
totalR = [];


for i = 1:L

    reflection_coefficients = 1*exp(1j*transpose((b-a).*rand(M,1)+a)*(2*pi)/360); %% Reflection coefficients vector
    

    A = zeros(Nr,Nt);

    for l = 1:M
        A1 = reflection_coefficients(l)*R(:,l)*T(l,:);
        A = A + A1; 
    end

   
    %%% Generation of the channel with IRS

    H_tilde = H + A;
    D2 = rank(H_tilde);
    
    %%% Singular value decomposition and acquisition of singular values

    [U,D,V] = svd(H_tilde); 
    lam = diag(D).^2; %Singular Values 


    %%% Water-filling power allocation

    wline = wfill(1./lam, P_max_linear,accuracy);
    Pi = max(wline-1./lam,0); % Tx power allocated to each Tx antenna
    power = zeros(1,size(V,2));
    power(1:length(lam)) = Pi; % optimal Tx power per Tx antenna + added zero to fit dimension of V;

    %%% Generation of the transmit covariance matrix

    Q = V*diag(power)*V'; % Tx covariance matrix

    %%% Channel capacity with IRS 

    Capacity_IRS = log2(det(eye(Nr)+1/sigma_lineal*H_tilde*Q*H_tilde'));
    
    totalQ = [totalQ;Q];
    totalC = [totalC Capacity_IRS];
    totalR = [totalR; reflection_coefficients];
    
end

[optCvalue, optCindex] = max(totalC);

%%% Selecting the optimal Q and reflection coefficients giving the largest
%%% rate

optQ = [totalQ(optCindex*4-3,:);totalQ(optCindex*4-2,:);totalQ(optCindex*4-1,:);totalQ(optCindex*4,:)];
optR = totalR(optCindex,:);
optCapacity = optCvalue;

end



