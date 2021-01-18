function [Sel_Mat_2D, K_miu_MS_Re, K_miu_MS_Im, K_niu_MS_Re, K_niu_MS_Im] = Spatial_Smoothing_2D(N_sub, G_h, G_v)

    N_h_sub = N_sub(1);
    N_v_sub = N_sub(2);
    M_h_sub = N_h_sub - G_h + 1;
    M_v_sub = N_v_sub - G_v + 1;
    M_sub = M_h_sub*M_v_sub;
    M = N_h_sub*N_v_sub;
    G_g = G_h*G_v;
    Sel_Mat_2D = zeros(M_sub, M, G_g);
    for g_h = 1:G_h
        Sel_Mat_h = [zeros(M_h_sub,g_h-1), eye(M_h_sub), zeros(M_h_sub,G_h-g_h)];
        for g_v = 1:G_v
            Sel_Mat_v = [zeros(M_v_sub,g_v-1), eye(M_v_sub), zeros(M_v_sub,G_v-g_v)];
            g_g = (g_h-1)*G_v + g_v;
            Sel_Mat_2D(:,:,g_g) = kron(Sel_Mat_v,Sel_Mat_h);
        end
    end
    J_M_miu_BS_2 = Select_M2(M_h_sub);
    J_M_niu_BS_2 = Select_M2(M_h_sub);
    J_miu_BS_2 = kron(eye(M_h_sub),J_M_miu_BS_2);
    J_niu_BS_2 = kron(J_M_niu_BS_2,eye(M_h_sub));
    K_miu_BS = Q(size(J_miu_BS_2,1))'*J_miu_BS_2*Q(size(J_miu_BS_2,2));
    K_niu_BS = Q(size(J_niu_BS_2,1))'*J_niu_BS_2*Q(size(J_niu_BS_2,2));
    K_miu_MS_Re = 2*real(K_miu_BS); K_miu_MS_Im = 2*imag(K_miu_BS);
    K_niu_MS_Re = 2*real(K_niu_BS); K_niu_MS_Im = 2*imag(K_niu_BS);

end