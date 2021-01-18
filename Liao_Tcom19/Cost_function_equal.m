function Phi_com_new = Cost_function_equal(theta_angle, Phi_com, N_miu, p, q, Lp)

    n_angle = length(theta_angle);
    if n_angle == 1
        JM_pq = Jacobi_Matrix(Lp, p, q, theta_angle);
    else
        JM_pq = zeros(Lp,Lp,n_angle);
        for i_angle = 1:n_angle
            angle = theta_angle(i_angle);
            JM_pq(:,:,i_angle) = Jacobi_Matrix(Lp, p, q, angle);
        end
    end
    tri_dem = size(JM_pq,3);
    Delta_psi_JM_pq_temp = zeros(tri_dem,N_miu);
    Phi_n_dot_com = zeros(Lp,Lp,N_miu,tri_dem);
    for i_dem = 1:tri_dem
        JM_pq_sel = JM_pq(:,:,i_dem);
        for n_miu = 1:N_miu
            Phi_n = Phi_com(:,:,n_miu);
            Phi_n_dot = JM_pq_sel.'*Phi_n*JM_pq_sel;
            Delta_psi_JM_pq_temp(i_dem,n_miu) = norm(tril(Phi_n_dot,-1),'fro')^2 - norm(tril(Phi_n,-1),'fro')^2;
            Phi_n_dot_com(:,:,n_miu,i_dem) = Phi_n_dot;
        end
    end
    Delta_psi_JM_pq = sum(Delta_psi_JM_pq_temp,2);
    [min_val, index] = min(Delta_psi_JM_pq);
    if min_val < 0
        Phi_com_new = Phi_n_dot_com(:,:,:,index);
    else
        Phi_com_new = Phi_com;
    end

end

function Jac_Mat = Jacobi_Matrix(n, p, q, angle)
    In = eye(n);
    c = cos(angle);
    s = sin(angle);
    if q > p
        In(p,p) = c;
        In(q,q) = c;
        In(p,q) = s;
        In(q,p) = -s;
        Jac_Mat = In;
    else
        error('Error: q should be > p!');
    end
end
