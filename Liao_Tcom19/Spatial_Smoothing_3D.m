function [Sel_Mat_3D, K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, K_tau_Re, K_tau_Im] = Spatial_Smoothing_3D(N_sub, K, G_h, G_v, G_k)

    N_h_sub = N_sub(1);
    N_v_sub = N_sub(2);
    M_h_sub = N_h_sub - G_h + 1;
    M_v_sub = N_v_sub - G_v + 1;
    K_sub = K - G_k +1;
    M_sub = M_h_sub*M_v_sub*K_sub;
    M = N_h_sub*N_v_sub*K;
    G_g = G_h*G_v*G_k;
    Sel_Mat_3D = zeros(M_sub, M, G_g);
    for g_h = 1:G_h
        Sel_Mat_h = [zeros(M_h_sub,g_h-1), eye(M_h_sub), zeros(M_h_sub,G_h-g_h)];
        for g_v = 1:G_v
            Sel_Mat_v = [zeros(M_v_sub,g_v-1), eye(M_v_sub), zeros(M_v_sub,G_v-g_v)];
            for g_k = 1:G_k
                Sel_Mat_k = [zeros(K_sub,g_k-1), eye(K_sub), zeros(K_sub,G_k-g_k)];
                g_g = (g_h-1)*G_v*G_k + (g_v-1)*G_k + g_k;
                Sel_Mat_3D(:,:,g_g) = kron(Sel_Mat_k,kron(Sel_Mat_v,Sel_Mat_h));
            end
        end
    end
    J_M_miu_BS_2 = Select_M2(M_h_sub);
    J_M_niu_BS_2 = Select_M2(M_h_sub);
    J_M_tau_2 = Select_M2(K_sub);
    J_miu_BS_2 = kron(eye(K_sub*M_h_sub),J_M_miu_BS_2);
    J_niu_BS_2 = kron(kron(eye(K_sub),J_M_niu_BS_2),eye(M_h_sub));
    J_tau_2 = kron(J_M_tau_2,eye(M_h_sub*M_h_sub));
    K_miu_BS = Q(size(J_miu_BS_2,1))'*J_miu_BS_2*Q(size(J_miu_BS_2,2));
    K_niu_BS = Q(size(J_niu_BS_2,1))'*J_niu_BS_2*Q(size(J_niu_BS_2,2));
    K_tau = Q(size(J_tau_2,1))'*J_tau_2*Q(size(J_tau_2,2));
    K_miu_BS_Re = 2*real(K_miu_BS); K_miu_BS_Im = 2*imag(K_miu_BS);
    K_niu_BS_Re = 2*real(K_niu_BS); K_niu_BS_Im = 2*imag(K_niu_BS);
    K_tau_Re = 2*real(K_tau); K_tau_Im = 2*imag(K_tau);

end