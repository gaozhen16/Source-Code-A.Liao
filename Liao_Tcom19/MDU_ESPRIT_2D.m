function N_miu_est = MDU_ESPRIT_2D(Y_d, Lp, K_miu_Re, K_miu_Im, K_niu_Re, K_niu_Im, Sel_Mat_2D, N_sub, G_h, G_v)

    N_h_sub = N_sub(1);
    N_v_sub = N_sub(2);
    M_h_sub = N_h_sub - G_h + 1;
    M_v_sub = N_v_sub - G_v + 1;
    M_sub = M_h_sub*M_v_sub;
    G_g = G_h*G_v;
    N_col = size(Y_d,2);
    N_sub = N_col*G_g;
    Y_d_smoothed = zeros(M_sub, N_sub);
    for gg = 1:G_g
        Y_d_smoothed(:,(gg-1)*N_col+1:gg*N_col) = Sel_Mat_2D(:,:,gg)*Y_d;
    end
    [M, N] = size(Y_d_smoothed);
    Y_d_real_temp = [Y_d_smoothed, flipud(eye(M))*conj(Y_d_smoothed)*flipud(eye(N))];
    Y_d_real = real(Q(M)'*Y_d_real_temp*Q(2*N));
    [E_s, ~, ~] = svds(Y_d_real, Lp);
    U_miu_1 = K_miu_Re*E_s; U_miu_2 = K_miu_Im*E_s; E_miu = [U_miu_1, U_miu_2];
    U_niu_1 = K_niu_Re*E_s; U_niu_2 = K_niu_Im*E_s; E_niu = [U_niu_1, U_niu_2];
    [U_miu, ~] = eig(E_miu'*E_miu); E_miu_12 = U_miu(1:Lp,1:Lp); E_miu_22 = U_miu(Lp+1:end,1:Lp);
    [U_niu, ~] = eig(E_niu'*E_niu); E_niu_12 = U_niu(1:Lp,1:Lp); E_niu_22 = U_niu(Lp+1:end,1:Lp);
    Phi_miu = - E_miu_12/E_miu_22;
    Phi_niu = - E_niu_12/E_niu_22;
    Psi_Joint = Phi_miu +1j*Phi_niu;

    Lambda = eig(Psi_Joint);
    miu_est = 2*atan(real(Lambda));
    niu_est = 2*atan(imag(Lambda));
    N_miu_est = [miu_est, niu_est];

end
