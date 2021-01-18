function N_miu_est = Ad_Hoc_Approach(Phi_com, Lp)
% An ad hoc approach for pairing with different paths
    N_miu = size(Phi_com,3);
    N_miu_est = zeros(Lp,N_miu);
    order_N_miu = N_miu:-1:1;
    for i_miu = 1:N_miu
        n_miu = order_N_miu(i_miu);
        if n_miu == N_miu
            [T_N_miu, Omega_N_tilde] = eig(Phi_com(:,:,n_miu));
            N_miu_est(:,n_miu) = 2*atan(diag(Omega_N_tilde));
        else
            Omega_n_tilde = T_N_miu\Phi_com(:,:,n_miu)*T_N_miu;
            N_miu_est(:,n_miu) = 2*atan(diag(Omega_n_tilde));
        end
    end
end