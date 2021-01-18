function N_miu_est = SSD_Algorithm(Phi_com, N_sweeps, Lp)

    N_miu = size(Phi_com,3);

    if Lp > 1
    Combination_temp = [];
    for i_q = 2:Lp
        for i_p = 1:i_q-1
            Combination_temp = [Combination_temp,i_q, i_p];
        end
    end
    N_com = nchoosek(Lp,2);
    Combination = zeros(N_com,2);
    for n_com1 = 1:N_com
        Combination(n_com1,:) = Combination_temp((n_com1-1)*2+1:n_com1*2);
    end

    for n_sweeps = 1:N_sweeps
        for n_com2 = 1:N_com
            n_pq = Combination(n_com2,:);
            p = n_pq(2);
            q = n_pq(1);
            theta_angle = Solve_theta_angle(Phi_com, N_miu, p, q);
            if theta_angle == 0
                Phi_com_new = Phi_com;
            else
                Phi_com_new = cost_function_equal(theta_angle, Phi_com, N_miu, p, q, Lp);
            end
            Phi_com = Phi_com_new;
        end
    end

    end

    U_com_last = Phi_com;
    N_miu_est = zeros(Lp,N_miu);
    for n_miu = 1:N_miu
        U_n_miu = U_com_last(:,:,n_miu);
        N_miu_est(:,n_miu) = 2*atan(diag(U_n_miu));
    end

end
