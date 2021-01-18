function N_miu_est = MDU_ESPRIT_3D(Y_u, Lp, K_miu_Re, K_miu_Im, K_niu_Re, K_niu_Im, K_tau_Re, K_tau_Im, Sel_Mat_3D, N_sub, K, G_h, G_v, G_k, method)

    N_h_sub = N_sub(1);
    N_v_sub = N_sub(2);
    M_h_sub = N_h_sub - G_h + 1;
    M_v_sub = N_v_sub - G_v + 1;
    K_sub = K - G_k +1;
    M_sub = M_h_sub*M_v_sub*K_sub;
    G_g = G_h*G_v*G_k;
    N_col = size(Y_u,2);
    N_sub = N_col*G_g;
    Y_u_smoothed = zeros(M_sub, N_sub);
    for gg = 1:G_g
        Y_u_smoothed(:,(gg-1)*N_col+1:gg*N_col) = Sel_Mat_3D(:,:,gg)*Y_u;
    end
    [M, N] = size(Y_u_smoothed);
    Y_u_real_temp = [Y_u_smoothed, flipud(eye(M))*conj(Y_u_smoothed)*flipud(eye(N))];
    Y_u_real = real(Q(M)'*Y_u_real_temp*Q(2*N));
    [E_s, ~, ~] = svds(Y_u_real, Lp);
    U_miu_1 = K_miu_Re*E_s; U_miu_2 = K_miu_Im*E_s; E_miu = [U_miu_1, U_miu_2];
    U_niu_1 = K_niu_Re*E_s; U_niu_2 = K_niu_Im*E_s; E_niu = [U_niu_1, U_niu_2];
    U_tau_1 = K_tau_Re*E_s; U_tau_2 = K_tau_Im*E_s; E_tau = [U_tau_1, U_tau_2];
    [U_miu, ~] = eig(E_miu'*E_miu); E_miu_12 = U_miu(1:Lp,1:Lp); E_miu_22 = U_miu(Lp+1:end,1:Lp);
    [U_niu, ~] = eig(E_niu'*E_niu); E_niu_12 = U_niu(1:Lp,1:Lp); E_niu_22 = U_niu(Lp+1:end,1:Lp);
    [U_tau, ~] = eig(E_tau'*E_tau); E_tau_12 = U_tau(1:Lp,1:Lp); E_tau_22 = U_tau(Lp+1:end,1:Lp);
    Phi_miu = - E_miu_12/E_miu_22;
    Phi_niu = - E_niu_12/E_niu_22;
    Phi_tau = - E_tau_12/E_tau_22;
    N_miu = 3;
    Phi_com = zeros(Lp,Lp,N_miu);
    Phi_com(:,:,1) = Phi_miu;
    Phi_com(:,:,2) = Phi_niu;
    Phi_com(:,:,3) = Phi_tau;
    if method(1) == 'A'
        N_miu_est = Ad_Hoc_Approach(Phi_com, Lp);
    else
        N_sweeps = 10;
        N_miu_est = SSD_Algorithm(Phi_com, N_sweeps, Lp);
    end

end