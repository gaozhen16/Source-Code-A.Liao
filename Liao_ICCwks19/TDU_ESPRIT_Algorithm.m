function N_miu_est = TDU_ESPRIT_Algorithm(Y_bar, Lp, K_miu_BS_Re, K_niu_BS_Re, K_miu_tau_Re, K_miu_BS_Im, K_niu_BS_Im, K_miu_tau_Im, M_sub, G_g, Sel_Mat_3D)

N_column = size(Y_bar,2);
N_sub = N_column*G_g;
Y_bar_SS = zeros(M_sub, N_sub);
for g_g = 1:G_g
    Y_bar_SS(:,(g_g-1)*N_column+1:g_g*N_column) = Sel_Mat_3D(:,:,g_g)*Y_bar;
end
[M_ss, N_ss] = size(Y_bar_SS);
H_bar_real_temp = [Y_bar_SS, flipud(eye(M_ss))*conj(Y_bar_SS)*flipud(eye(N_ss))];
H_bar_real = real(Q(M_ss)'*H_bar_real_temp*Q(2*N_ss));
[E_s, ~, ~] = svds(H_bar_real, Lp);
U_miu_BS_1 = K_miu_BS_Re*E_s; U_miu_BS_2 = K_miu_BS_Im*E_s; E_miu_BS = [U_miu_BS_1, U_miu_BS_2];
U_niu_BS_1 = K_niu_BS_Re*E_s; U_niu_BS_2 = K_niu_BS_Im*E_s; E_niu_BS = [U_niu_BS_1, U_niu_BS_2];
U_miu_tau_1 = K_miu_tau_Re*E_s; U_miu_tau_2 = K_miu_tau_Im*E_s; E_miu_tau = [U_miu_tau_1, U_miu_tau_2];
[U_miu_BS, ~] = eig(E_miu_BS'*E_miu_BS); E_miu_BS_12 = U_miu_BS(1:Lp,1:Lp); E_miu_BS_22 = U_miu_BS(Lp+1:end,1:Lp);
[U_niu_BS, ~] = eig(E_niu_BS'*E_niu_BS); E_niu_BS_12 = U_niu_BS(1:Lp,1:Lp); E_niu_BS_22 = U_niu_BS(Lp+1:end,1:Lp);
[U_miu_tau, ~] = eig(E_miu_tau'*E_miu_tau); E_miu_tau_12 = U_miu_tau(1:Lp,1:Lp); E_miu_tau_22 = U_miu_tau(Lp+1:end,1:Lp);
Phi_miu_BS = - E_miu_BS_12/E_miu_BS_22;
Phi_niu_BS = - E_niu_BS_12/E_niu_BS_22;
Phi_miu_tau = - E_miu_tau_12/E_miu_tau_22;
N_miu = 3;
Phi_com = zeros(Lp,Lp,N_miu);
Phi_com(:,:,1) = Phi_miu_BS;
Phi_com(:,:,2) = Phi_niu_BS;
Phi_com(:,:,3) = Phi_miu_tau;
N_sweeps = 10;
N_miu_est = SSD_Algorithm(Phi_com, N_sweeps, Lp);

end
