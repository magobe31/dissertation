
     %%%% Genaration of the direct channel H between the transmitter and the
     %%%% receiver

function [H] = direct_channel(H_hat,dD_hat,Nr,Nt,lambda,theta_AoA,theta_AoD,Kd,do,beta_o_lineal,dA,alpha_d)

ar1 = [];
at1 = [];

for i = 1:Nr
    ar = exp(i*2*pi*(i-1)*dA*(sin(theta_AoA))/lambda);
    ar1 = [ar1 ar];
end

for i = 1:Nt
    at = exp(i*2*pi*(i-1)*dA*(sin(theta_AoD))/lambda);
    at1 = [at1 at];
end

ar1 = ar1';
at1 = at1';

H_LoS = ar1*at1';

H_NLoS_direct = sqrt(1/2)*(randn(Nr,Nt)+1j*randn(Nr,Nt)); %% Rayleigh channel matrix of the direct link between Tx-Rx (variance of 1)

dD = sqrt(dD_hat^2+H_hat^2);
beta_d = beta_o_lineal*(dD/do)^-alpha_d;

H = sqrt(beta_d/(Kd+1))*(sqrt(Kd)*H_LoS+H_NLoS_direct);

end
