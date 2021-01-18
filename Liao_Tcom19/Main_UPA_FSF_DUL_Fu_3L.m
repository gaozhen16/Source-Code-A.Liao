clc, clear, warning off;

N_ant = 12;
N_ms = [N_ant,N_ant];
N_bs = [N_ant,N_ant];
N_MS = N_ms(1)*N_ms(2);
N_BS = N_bs(1)*N_bs(2);
FFT_len = 128;
K = FFT_len;
BW = 200e6;
fs = BW;
CIR_max_L = 16;
fc = 30e9;
lambda = 3e8/fc;
d_ant = lambda/2;
sigma_2_alpha = 1;
awgn_en = 1;
method = 'SSD';
N_bits_pre = 3;
N_bits_fb = 10;
N_Bits_fb = 2^N_bits_fb;
Lp = 3;
N_OFDM = 3;
N_RF_BS = 4;
Ns_u = N_RF_BS-1;
N_bs_h_sub = 8;
N_bs_v_sub = N_bs_h_sub;
N_BS_sub = [N_bs_h_sub, N_bs_v_sub];
N_u_B = ceil(N_bs_h_sub*N_bs_v_sub/Ns_u);
N_bs_re = N_u_B*Ns_u - N_bs_h_sub*N_bs_v_sub;
Sel_u = [eye(N_bs_h_sub*N_bs_v_sub),zeros(N_bs_h_sub*N_bs_v_sub,N_bs_re)];
N_RF_MS = N_RF_BS;
Ns_d = N_RF_MS - 1;
N_ms_h_sub = N_bs_h_sub;
N_ms_v_sub = N_bs_h_sub;
N_MS_sub = [N_ms_h_sub, N_ms_v_sub];
N_d_B = ceil(N_ms_h_sub*N_ms_v_sub/Ns_d);
N_ms_re = N_d_B*Ns_d - N_ms_h_sub*N_ms_v_sub;
Sel_d = [eye(N_ms_h_sub*N_ms_v_sub),zeros(N_ms_h_sub*N_ms_v_sub,N_ms_re)];
T_pilot = N_OFDM*(N_u_B + N_d_B);
G_h_2 = 2; G_v_2 = 2;
[Sel_Mat_2D, K_miu_MS_Re, K_miu_MS_Im, K_niu_MS_Re, K_niu_MS_Im] = Spatial_Smoothing_2D(N_MS_sub, G_h_2, G_v_2);
G_h_3 = 2; G_v_3 = 2; G_k = K/2;
[Sel_Mat_3D, K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, K_tau_Re, K_tau_Im] = Spatial_Smoothing_3D(N_BS_sub, K, G_h_3, G_v_3, G_k);
[W_tilde_u, W_Mat_u] = Combiner_Design(N_bs, N_BS_sub, N_RF_MS, Ns_u, N_bits_pre);
W_H_blkdiag_u = zeros(N_u_B*Ns_u,N_u_B*N_BS);
for n_u = 1:N_u_B
    W_H_blkdiag_u((n_u-1)*Ns_u+1:n_u*Ns_u,(n_u-1)*N_BS+1:n_u*N_BS) = W_Mat_u(:,:,n_u)';
end
[W_tilde_d, W_Mat_d] = Combiner_Design(N_ms, N_MS_sub, N_RF_BS, Ns_d, N_bits_pre);
W_H_blkdiag_d = zeros(N_d_B*Ns_d,N_d_B*N_MS);
for n_d = 1:N_d_B
	W_H_blkdiag_d((n_d-1)*Ns_d+1:n_d*Ns_d,(n_d-1)*N_MS+1:n_d*N_MS) = W_Mat_d(:,:,n_d)';
end
iterMax = 1e3;
SNR_dBs = -15:5:10;
NMSE = zeros(1,length(SNR_dBs));

for ii = 1:length(SNR_dBs)   % 
    tic
    sigma2 = 10^(-(SNR_dBs(ii)/10));
    sigma = sqrt(sigma2);
    for iter = 1:iterMax
        [H_f, theta_BS, phi_BS, tau, theta_MS, phi_MS, miu_BS, niu_BS, miu_tau, miu_MS, niu_MS, alpha] = ...
            FSF_Downlink_Channel(N_ms, N_bs, fc, Lp, sigma_2_alpha, fs, K, CIR_max_L);
        S_d = exp(-1i*2*pi*rand(Ns_d, N_OFDM));
        F_RF_d = exp(1j*2*pi*rand(N_BS,N_RF_MS));
        Quantized_F_RF_d = Quantize(F_RF_d, N_bits_pre);
        F_BB_d = exp(1j*2*pi*rand(N_RF_MS,Ns_d));
        F_temp_d = Quantized_F_RF_d*F_BB_d;
        lambda_d = sqrt(N_RF_MS)/norm(F_temp_d,'fro');
        F_d = lambda_d*F_temp_d;
        Y_d = zeros(N_ms_h_sub*N_ms_v_sub,N_OFDM*K);
        Y_d_com = zeros(N_ms_h_sub*N_ms_v_sub, N_OFDM, K);
        for kk1 = 1:K
            Y_DL_k_temp = W_tilde_d'*H_f(:,:,kk1)*F_d*S_d + ...
                awgn_en*sigma*W_H_blkdiag_d*(normrnd(0,1,N_d_B*N_MS,N_OFDM) + 1i*normrnd(0,1,N_d_B*N_MS,N_OFDM))/sqrt(2);
            Y_DL_k = Sel_d*Y_DL_k_temp;
            Y_d(:,(kk1-1)*N_OFDM+1:kk1*N_OFDM) = Y_DL_k;
            Y_d_com(:,:,kk1) = Y_DL_k;
        end
        Lp_est = Lp;
        N_miu_est_d = MDU_ESPRIT_2D(Y_d, Lp_est, K_miu_MS_Re, K_miu_MS_Im, K_niu_MS_Re, K_niu_MS_Im, Sel_Mat_2D, N_MS_sub, G_h_2, G_v_2);
        miu_MS_est = N_miu_est_d(:,1);
       	niu_MS_est = N_miu_est_d(:,2);
        m_MS = (0:N_ms(1)-1).';
       	n_MS = (0:N_ms(2)-1).';
       	A_miu_MS_est = exp(1i*m_MS*miu_MS_est.')/sqrt(N_ms(1));
       	A_niu_MS_est = exp(1i*n_MS*niu_MS_est.')/sqrt(N_ms(2));
       	A_MS_est = Khatri_Rao(A_niu_MS_est,A_miu_MS_est);
       	phi_MS_hat = asin(niu_MS_est*lambda/(2*pi*d_ant));
       	theta_MS_hat = asin(miu_MS_est*lambda./(2*pi*d_ant*cos(phi_MS_hat)));
       	phi_MS_hat_quan = round(phi_MS_hat*N_Bits_fb/(2*pi/3)) *(2*pi/3)/N_Bits_fb;
      	theta_MS_hat_quan = round(theta_MS_hat*N_Bits_fb/(2*pi/3)) *(2*pi/3)/N_Bits_fb;
       	miu_MS_fb = 2*pi/lambda*d_ant*(sin(theta_MS_hat_quan).*cos(phi_MS_hat_quan));
       	niu_MS_fb = 2*pi/lambda*d_ant*sin(phi_MS_hat_quan);
        if Lp_est > N_RF_MS
            Num_ant = floor(N_RF_MS*N_MS/Lp_est);
            Num_ant_re = mod(N_RF_MS*N_MS,Lp_est);
            Num_ant_set = Num_ant*ones(Lp_est,1);
            if Num_ant_re ~= 0
                Num_ant_set(1:Num_ant_re) = Num_ant_set(1:Num_ant_re) + ones(Num_ant_re,1);
            end
            Index_Sum = (1:N_RF_MS*N_MS).';
            F_RF_u_vec = zeros(N_RF_MS*N_MS,1);
            for llp = 1:Lp_est
                Index_temp1 = Index_Sum(sum(Num_ant_set(1:(llp-1)))+1:sum(Num_ant_set(1:llp)));
                Index_temp2 = mod(Index_temp1,N_MS);
                Index_temp2(Index_temp2 == 0) = N_MS;
                F_RF_u_vec(Index_temp1) = A_MS_est(Index_temp2,llp);
            end
            F_RF_u = reshape(F_RF_u_vec,N_MS,N_RF_MS);
            Quantized_F_RF_u = Quantize(conj(F_RF_u), N_bits_pre);
            F_BB_u = exp(1j*2*pi*rand(N_RF_MS,Ns_d));
            F_temp_u = Quantized_F_RF_u*F_BB_u;
            lambda_u = sqrt(N_RF_MS)/norm(F_temp_u,'fro');
            F_u = lambda_u*F_temp_u;
        else
            if mod(N_RF_MS,Lp_est) == 0
                N_repmat = N_RF_MS/Lp_est;
                F_RF_u = repmat(A_MS_est,1,N_repmat);
            else
                N_repmat = floor(N_RF_MS/Lp_est);
                F_RF_u = zeros(N_MS,N_RF_MS);
                F_RF_u(:,1:Lp_est*N_repmat) = repmat(A_MS_est,1,N_repmat);
                N_RF_u_re = mod(N_RF_MS,Lp_est);
                Num_ant = floor(N_RF_u_re*N_MS/Lp_est);
                Num_ant_re = mod(N_RF_u_re*N_MS,Lp_est);
                Num_ant_set = Num_ant*ones(Lp_est,1);
                if Num_ant_re ~= 0
                    Num_ant_set(1:Num_ant_re) = Num_ant_set(1:Num_ant_re) + ones(Num_ant_re,1);
                end
                Index_Sum = (1:N_RF_u_re*N_MS).';
                F_RF_u_vec = zeros(N_RF_u_re*N_MS,1);
                for llp = 1:Lp_est
                    Index_temp1 = Index_Sum(sum(Num_ant_set(1:(llp-1)))+1:sum(Num_ant_set(1:llp)));
                    Index_temp2 = mod(Index_temp1,N_MS);
                    Index_temp2(Index_temp2 == 0) = N_MS;
                    F_RF_u_vec(Index_temp1) = A_MS_est(Index_temp2,llp);
                end
                F_RF_u(:,Lp_est*N_repmat+1:end) = reshape(F_RF_u_vec,N_MS,N_RF_u_re);
            end
            Quantized_F_RF_u = Quantize(conj(F_RF_u), N_bits_pre);
            F_BB_u = exp(1j*2*pi*rand(N_RF_MS,Ns_d));
            F_temp_u = Quantized_F_RF_u*F_BB_u;
            lambda_u = sqrt(N_RF_MS)/norm(F_temp_u,'fro');
            F_u = lambda_u*F_temp_u;
        end
        S_u = exp(-1i*2*pi*rand(Ns_u, N_OFDM));
        Y_u = zeros(N_bs_h_sub*N_bs_v_sub*N_OFDM,K);
        for kk2 = 1:K
            Y_u_k_temp = W_tilde_u'*H_f(:,:,kk2).'*F_u*S_u + ...
                awgn_en*sigma*W_H_blkdiag_u*(normrnd(0,1,N_u_B*N_BS,N_OFDM) + 1i*normrnd(0,1,N_u_B*N_BS,N_OFDM))/sqrt(2);
            Y_u_k = Sel_u*Y_u_k_temp;
            Y_u(:,kk2) = reshape(Y_u_k.',N_bs_h_sub*N_bs_v_sub*N_OFDM,1);
        end
        y_u_tilde = reshape(Y_u,K*N_bs_h_sub*N_bs_v_sub*N_OFDM,1);
        Y_u_bar = reshape(y_u_tilde,N_OFDM,K*N_bs_h_sub*N_bs_v_sub).';
        N_miu_u_est = MDU_ESPRIT_3D(Y_u_bar, Lp_est, K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, K_tau_Re, K_tau_Im,...
            Sel_Mat_3D, N_BS_sub, K, G_h_3, G_v_3, G_k, method);
        if Lp_est == 1
            miu_BS_est = - N_miu_u_est(:,1);
        	niu_BS_est = - N_miu_u_est(:,2);
            tau_est = N_miu_u_est(:,3);
        else
            [tau_est, Index] = sort(N_miu_u_est(:,3),'descend');
            niu_BS_est_temp = - N_miu_u_est(:,2);
            miu_BS_est_temp = - N_miu_u_est(:,1);
            niu_BS_est = niu_BS_est_temp(Index);
            miu_BS_est = miu_BS_est_temp(Index);
        end
        k_tau = (0:(K-1))';
        A_tau_est = exp(1i*k_tau*tau_est.');
        m_BS = (0:N_bs(1)-1).';
        n_BS = (0:N_bs(2)-1).';
        A_miu_BS_est = exp(1i*m_BS*miu_BS_est.')/sqrt(N_bs(1));
        A_niu_BS_est = exp(1i*n_BS*niu_BS_est.')/sqrt(N_bs(2));
        A_BS_est = Khatri_Rao(A_niu_BS_est,A_miu_BS_est);
        if Lp_est == 1
            m_MS = (0:N_ms(1)-1).';
            n_MS = (0:N_ms(2)-1).';
            A_miu_MS_fb = exp(1i*m_MS*miu_MS_fb.')/sqrt(N_ms(1));
            A_niu_MS_fb = exp(1i*n_MS*niu_MS_fb.')/sqrt(N_ms(2));
            A_MS_fb = Khatri_Rao(A_niu_MS_fb,A_miu_MS_fb);
            A_temp = Khatri_Rao(Khatri_Rao(A_tau_est,(Sel_u*W_tilde_u'*conj(A_BS_est))),(A_MS_fb.'*F_u*S_u).');
            alpha_est_temp = A_temp\y_u_tilde;
            alpha_est = alpha_est_temp/sqrt(N_MS*N_BS/Lp_est);
        else
            Permutation = perms(Lp_est:-1:1); N_perms = size(Permutation,1);
            alpha_est_Mat = zeros(Lp_est,N_perms);
            Difference = zeros(N_perms,1);
            for n_p = 1:N_perms
                index_per = Permutation(n_p,:).';
                miu_MS_est_temp = miu_MS_fb(index_per);
                niu_MS_est_temp = niu_MS_fb(index_per);
                A_miu_MS_fb_temp = exp(1i*m_MS*miu_MS_est_temp.')/sqrt(N_ms(1));
                A_niu_MS_fb_temp = exp(1i*n_MS*niu_MS_est_temp.')/sqrt(N_ms(2));
                A_MS_fb_temp = Khatri_Rao(A_niu_MS_fb_temp,A_miu_MS_fb_temp);
                A_temp = Khatri_Rao(Khatri_Rao(A_tau_est,(Sel_u*W_tilde_u'*conj(A_BS_est))),(A_MS_fb_temp.'*F_u*S_u).');
                alpha_est_temp = A_temp\y_u_tilde;
                y_u_est = A_temp*alpha_est_temp;
                alpha_est_Mat(:,n_p) = alpha_est_temp;
                Difference(n_p) = norm(y_u_tilde - y_u_est)^2;
            end
            [~,N_p] = min(Difference);
            index_order = Permutation(N_p,:).';
            miu_MS_est = miu_MS_fb(index_order);
            niu_MS_est = niu_MS_fb(index_order);
            alpha_est = alpha_est_Mat(:,N_p)/sqrt(N_MS*N_BS/Lp_est);
            A_miu_MS_est = exp(1i*m_MS*miu_MS_est.')/sqrt(N_ms(1));
            A_niu_MS_est = exp(1i*n_MS*niu_MS_est.')/sqrt(N_ms(2));
            A_MS_est = Khatri_Rao(A_niu_MS_est,A_miu_MS_est);
        end
        H_f_est = zeros(N_MS,N_BS,K);
        for kk_1 = 1:K
            D_diag_est = sqrt(N_MS*N_BS/Lp_est)*diag(alpha_est.*exp(1i*(kk_1-1)*tau_est));
            H_f_est(:,:,kk_1) = A_MS_est*D_diag_est*A_BS_est';
        end
        difference_channel_MSE = zeros(K,1);
        true_channel_MSE = zeros(K,1);
        for kk_2 = 1:K
            difference_channel_MSE(kk_2) = norm(H_f_est(:,:,kk_2) - H_f(:,:,kk_2),'fro')^2;
            true_channel_MSE(kk_2) = norm(H_f(:,:,kk_2),'fro')^2;
        end
        NMSE_temp = sum(difference_channel_MSE)/sum(true_channel_MSE);
        NMSE(ii) = NMSE(ii) + NMSE_temp;
        if mod(iter,100) == 0
            toc
            disp(['Pilot = ' num2str(T_pilot) ', SNR = ' num2str(SNR_dBs(ii)) ', iter = ' num2str(iter)...
                ' , NMSE_current = ' num2str(10*log10(NMSE_temp)) ' , NMSE_mean = ' num2str(10*log10(NMSE(ii)/iter))])
        end
    end
    toc
    disp(['Finished ',num2str(ii),'/', num2str(length(SNR_dBs)) ' , NMSE = ' num2str(10*log10(NMSE(ii)/iter))]);
end
NMSE = NMSE/iterMax;
NMSE_dB_ESPRIT_3L_Lp = 10*log10(NMSE);

disp('Finished all');

%% Plot
MarkerSize = 6;
LineWidth = 2;
Fontsize = 15;

figure
plot(SNR_dBs,NMSE_dB_ESPRIT_3L_Lp,'-.ko','MarkerFaceColor',[0 0 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('NMSE [dB]','Fontsize',Fontsize)
set(gca,'FontSize',Fontsize, 'linewidth',1.5);
set(gca, 'GridLineStyle', '-.');grid on;
set(gcf, 'position', [700 300 650 550]);
h1 = legend('{\it{L}}=3','Location','southwest');
set(h1,'Fontsize',14);
