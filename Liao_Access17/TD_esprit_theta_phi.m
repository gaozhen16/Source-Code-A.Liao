function [theta_est, phi_est] = TD_esprit_theta_phi(H_up_hat, m1, m2, Lp, lambda, d_ant)

    [M, N] = size(H_up_hat);
    if m1*M <= Lp
        error('m1 and N_BS is not large enough to separate all sources');
    end
    if m1 >= N
        error('m1 is too large');
    end
    if m2 >= M
        error('m2 is too large');
    end
    if ((Lp > (m1-1)*(M-m2+1)) || (Lp > m1*(M-m2)) || (Lp > 2*m2*(N-m1+1)))
        error('m1 or m2 need to change');
    end
    H_row = [];
    for i = 1:m2
       H_row = [H_row,H_up_hat(i:M-m2+i,:)];
    end
    M_m2 = M-m2+1;
    H_hankel = [];
    for j = 1:m1
    	H_col = [];
        for i = 1:m2;
            H_col = [H_col,H_row(:,(i-1)*N+j:i*N-m1+j)];
        end
        H_hankel = [H_hankel; H_col];
    end
    [H_h_row, H_h_col] = size(H_hankel);
    H_e = [H_hankel, flipud(eye(H_h_row))*conj(H_hankel)*flipud(eye(H_h_col))];
    T_L = kron(Q(m1)',Q(M_m2)');
    He_R = real(T_L*H_e*Q(2*H_h_col));
    [U_He_R, S_He_R, V_He_R] = svd(He_R);
    U_R = U_He_R(:,1:Lp);
    J2_phi  = [zeros(m1-1,1),eye(m1-1)];
    J2_theta  = [zeros(M_m2-1,1),eye(M_m2-1)];
    B_phi = kron(Q(m1-1)'*J2_phi*Q(m1),eye(M_m2));
    B_theta = kron(eye(m1),Q(M_m2-1)'*J2_theta*Q(M_m2));
    B_phi_R = real(B_phi);
    B_phi_I = imag(B_phi);
    B_theta_R = real(B_theta);
    B_theta_I = imag(B_theta);
    Ex_phi = B_phi_R*U_R;
    Ey_phi = B_phi_I*U_R;
    Ex_theta = B_theta_R*U_R;
    Ey_theta = B_theta_I*U_R;
    Joint_est = (Ex_phi)\Ey_phi+1i*(Ex_theta)\Ey_theta;
    Lambda = eig(Joint_est);
    Phi_tilde = real(Lambda);
    Theta_tilde = imag(Lambda);
    phi_est = -asin(atan(Phi_tilde)*lambda/(pi*d_ant));
    theta_est = -asin(atan(Theta_tilde)*lambda/(pi*d_ant));
end

function Q_N = Q(N)
if mod(N,2) == 0
    k = N/2;
    I_k = eye(k);
    J_k = flipud(eye(k));
    Q_N = [I_k,1j*I_k;...
           J_k,-1j*J_k]/sqrt(2);
else
    k = (N-1)/2;
    I_k = eye(k);
    J_k = flipud(eye(k));
    zero_col = zeros(k,1);
    Q_N = [I_k,zero_col,1j*I_k;...
           zero_col.',sqrt(2),zero_col.';...
           J_k,zero_col,-1j*J_k]/sqrt(2);
end
end