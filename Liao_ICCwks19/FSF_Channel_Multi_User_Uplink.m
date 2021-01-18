function [H_f, theta_BS, phi_BS, tau, miu_BS, niu_BS, miu_tau, alpha] = FSF_Channel_Multi_User_Uplink(N_bs, fc, N_U, sigma_2, fs, K, tao_max)

lambda = 3e8/fc;
d_ant = lambda/2;
az_ang_min = -pi/3;
az_ang_max = pi/3;
el_ang_min = -pi/3;
el_ang_max = pi/3;

%%
m_BS = (0:N_bs(1)-1).';
n_BS = (0:N_bs(2)-1).';

theta_BS = az_ang_min + (az_ang_max-az_ang_min).*rand(1,N_U);
phi_BS = el_ang_min + (el_ang_max-el_ang_min).*rand(1,N_U);
miu_BS = 2*pi/lambda*d_ant*(sin(theta_BS).*cos(phi_BS));
niu_BS = 2*pi/lambda*d_ant*sin(phi_BS);
A_miu_BS = exp(1i*m_BS*miu_BS)/sqrt(N_bs(1));
A_niu_BS = exp(1i*n_BS*niu_BS)/sqrt(N_bs(2));
A_BS = Khatri_Rao(A_niu_BS,A_miu_BS);

%%
tau_max = tao_max/fs;
tau = sort(tau_max.*rand(1,N_U));
miu_tau = -2*pi*fs*tau/K;

alpha_temp = sqrt(sigma_2/2)*(randn(1,N_U)+1i*randn(1,N_U));
alpha = sort(alpha_temp,'descend');

%%
N_BS = N_bs(1)*N_bs(2);
H_f = zeros(N_BS,N_U,K);
for kk = 1:K
  	H_f(:,:,kk) = A_BS*sqrt(N_BS/N_U)*diag(alpha.*exp(1i*(kk-1)*miu_tau));
end

end