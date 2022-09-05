
%%%% Generation of the Transmitter - IRS channel %%%%

function [T] = Transmitter_IRS_channel(Nt,M,dp_hat,dh_hat,dD_hat,do,alpha_TI,beta_o_lineal)
dTI = sqrt((dD_hat-dh_hat)^2+dp_hat^2);
beta_TI = beta_o_lineal*(dTI/do)^-alpha_TI;
NLoS_Tx_IRS = sqrt(1/2)*(randn(M,Nt)+1j*randn(M,Nt));
T = sqrt(beta_TI)*NLoS_Tx_IRS;
end


