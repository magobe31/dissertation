
%%%% Generation of the IRS - Receiver channel %%%%

function [R] = IRS_Receiver_channel(H_hat,Nr,dp_hat,dh_hat,do,alpha_IR,M,beta_o_lineal)

dIR = sqrt(dh_hat^2+dp_hat^2+H_hat^2);
beta_IR = beta_o_lineal*(dIR/do)^-alpha_IR;
NLoS_IRS_Rx = sqrt(1/2)*(randn(Nr,M)+1j*randn(Nr,M));
R = sqrt(beta_IR)*NLoS_IRS_Rx;
end

