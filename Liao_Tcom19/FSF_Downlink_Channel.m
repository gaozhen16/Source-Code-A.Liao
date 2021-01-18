function [H_f, theta_BS, phi_BS, tau, theta_MS, phi_MS, miu_BS, niu_BS, miu_tau, miu_MS, niu_MS, alpha] = ...
    FSF_Downlink_Channel(N_ms, N_bs, fc, Lp, sigma_2_alpha, fs, K, CIR_max_L)

    lambda = 3e8/fc;
    d_ant = lambda/2;
    az_ang_min = -pi/3;
    az_ang_max = pi/3;
    el_ang_min = -pi/3;
    el_ang_max = pi/3;

    check = length(N_ms)*length(N_bs);
    switch check
        case 1
            N_MS = N_ms;
            N_BS = N_bs;
            n_MS = (0:(N_MS-1)).';
            n_BS = (0:(N_BS-1)).';
            theta_MS = az_ang_min+(az_ang_max-az_ang_min).*rand(1,Lp);
            theta_BS = az_ang_min+(az_ang_max-az_ang_min).*rand(1,Lp);
            miu_MS = 2*pi/lambda*d_ant*sin(theta_MS);
            miu_BS = 2*pi/lambda*d_ant*sin(theta_BS);
            A_MS = exp(1i*n_MS*miu_MS)/sqrt(N_MS);
            A_BS = exp(1i*n_BS*miu_BS)/sqrt(N_BS);
            phi_BS = [];
            phi_MS = [];
            niu_BS = [];
            niu_MS = [];
        case 4
            m_MS = (0:N_ms(1)-1).';
            n_MS = (0:N_ms(2)-1).';
            N_MS = N_ms(1)*N_ms(2);
            theta_MS = az_ang_min + (az_ang_max-az_ang_min).*rand(1,Lp);
            phi_MS = el_ang_min + (el_ang_max-el_ang_min).*rand(1,Lp);
            miu_MS = 2*pi/lambda*d_ant*(sin(theta_MS).*cos(phi_MS));
            niu_MS = 2*pi/lambda*d_ant*sin(phi_MS);
            A_miu_MS = exp(1i*m_MS*miu_MS)/sqrt(N_ms(1));
            A_niu_MS = exp(1i*n_MS*niu_MS)/sqrt(N_ms(2));
            A_MS = Khatri_Rao(A_niu_MS,A_miu_MS);

            m_BS = (0:N_bs(1)-1).';
            n_BS = (0:N_bs(2)-1).';
            N_BS = N_bs(1)*N_bs(2);
            theta_BS = az_ang_min + (az_ang_max-az_ang_min).*rand(1,Lp);
            phi_BS = el_ang_min + (el_ang_max-el_ang_min).*rand(1,Lp);
            miu_BS = 2*pi/lambda*d_ant*(sin(theta_BS).*cos(phi_BS));
            niu_BS = 2*pi/lambda*d_ant*sin(phi_BS);
            A_miu_BS = exp(1i*m_BS*miu_BS)/sqrt(N_bs(1));
            A_niu_BS = exp(1i*n_BS*niu_BS)/sqrt(N_bs(2));
            A_BS = Khatri_Rao(A_niu_BS,A_miu_BS);
        otherwise
            error('Error: invalid anttena arguments!');
    end

    tau_max = CIR_max_L/fs;
    tau = sort(tau_max.*rand(1,Lp));
    miu_tau = -2*pi*fs*tau/K;
    alpha_temp = sqrt(sigma_2_alpha/2)*(randn(1,Lp)+1i*randn(1,Lp));
    alpha = sort(alpha_temp,'descend');

    H_f = zeros(N_MS,N_BS,K);
    for k = 1:K
        D_diag = sqrt(N_MS*N_BS/Lp)*diag(alpha.*exp(1i*(k-1)*miu_tau));
        H_f(:,:,k) = A_MS*D_diag*A_BS';
    end

end