clc, clear, warning off;
%-------------------------------------------------------------------------------------------------------------------------------------------------
% Multi-user FSF channel estimation in the uplink by expoliting TDU-ESPRIT algorithm to estimate the horizontal/vertical AoAs and delays at the BS.
%-------------------------------------------------------------------------------------------------------------------------------------------------

N_bs = [10,10];
N_BS = N_bs(1)*N_bs(2);
N_UE = 1;
N_RF = 4;
FFT_len = 128;
K = FFT_len;
BW = 200e6;
fs = BW;
tao_max = 16;
fc = 30e9;
lambda = 3e8/fc;
d_ant = lambda/2;
sigma_2_alpha = 1;
awgn_en = 1;
N_bits = 1;
N_U_set = 1:N_RF-1;
N_OFDM_set = [2,4,6];

iterMax = 2e3;
SNR_dBs = -10:5:15;
NMSE = zeros(length(N_U_set),length(SNR_dBs));

for uu = 1:length(N_U_set)
    N_U = N_U_set(uu);
    Ns = N_U;
    N_OFDM = N_OFDM_set(uu);
    N_h_sub = 6;
    N_v_sub = 6;
    N_hv_sub = N_h_sub*N_v_sub;
    N_BS_sub = [N_h_sub, N_v_sub];
    Nt = ceil(N_h_sub*N_v_sub/Ns);
    Nt_re = Nt*Ns - N_h_sub*N_v_sub;
    Sel_mat = [eye(N_h_sub*N_v_sub),zeros(N_h_sub*N_v_sub,Nt_re)];

    G_h = 2; G_v = 2; K_sub = K/2; G_k = K - K_sub +1;
    M_h_sub = N_h_sub - G_h + 1;
    M_v_sub = N_v_sub - G_v + 1;
    M_BS_sub = [M_h_sub,M_v_sub];
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

    J_M_miu_BS_2 = Select_M2(M_BS_sub(1));
    J_M_niu_BS_2 = Select_M2(M_BS_sub(2));
    J_M_miu_tau_2 = Select_M2(K_sub);
    J_miu_BS_2 = kron(eye(K_sub*M_BS_sub(2)),J_M_miu_BS_2);
    J_niu_BS_2 = kron(kron(eye(K_sub),J_M_niu_BS_2),eye(M_BS_sub(1)));
    J_miu_tau_2 = kron(J_M_miu_tau_2,eye(M_BS_sub(2)*M_BS_sub(1)));
    K_miu_BS = Q(size(J_miu_BS_2,1))'*J_miu_BS_2*Q(size(J_miu_BS_2,2));
    K_niu_BS = Q(size(J_niu_BS_2,1))'*J_niu_BS_2*Q(size(J_niu_BS_2,2));
    K_miu_tau = Q(size(J_miu_tau_2,1))'*J_miu_tau_2*Q(size(J_miu_tau_2,2));
    K_miu_BS_Re = 2*real(K_miu_BS); K_miu_BS_Im = 2*imag(K_miu_BS);
    K_niu_BS_Re = 2*real(K_niu_BS); K_niu_BS_Im = 2*imag(K_niu_BS);
    K_miu_tau_Re = 2*real(K_miu_tau); K_miu_tau_Im = 2*imag(K_miu_tau);

    [W_tilde_UL, W_Mat_UL] = Combiner_Design(N_bs, N_BS_sub, N_RF, Ns, N_bits);
    W_H_blkdiag_UL = zeros(Nt*Ns,Nt*N_BS);
    for n_t = 1:Nt
        W_H_blkdiag_UL((n_t-1)*Ns+1:n_t*Ns,(n_t-1)*N_BS+1:n_t*N_BS) = W_Mat_UL(:,:,n_t)';
    end

    for ii = 1:length(SNR_dBs)
        tic
        sigma2 = 10^(-(SNR_dBs(ii)/10));
        sigma = sqrt(sigma2);
        for iter = 1:iterMax
            [H_f, theta_BS, phi_BS, tau, miu_BS, niu_BS, miu_tau, alpha] = FSF_Channel_Multi_User_Uplink(N_bs, fc, N_U, sigma_2_alpha, fs, K, tao_max);
            S_pilot = exp(-1i*2*pi*rand(Ns, N_OFDM))/sqrt(Ns);
            Y_tilde = zeros(N_hv_sub*N_OFDM,K);
            for kk1 = 1:K
                Y_tilde_k = zeros(Nt*Ns,N_OFDM);
                for nnt = 1:Nt
                    Y_tilde_k((nnt-1)*Ns+1:nnt*Ns,:) = W_Mat_UL(:,:,nnt)'*(H_f(:,:,kk1)*S_pilot + ...
                        awgn_en*sigma*(normrnd(0,1,N_BS,N_OFDM) + 1i*normrnd(0,1,N_BS,N_OFDM))/sqrt(2));
                end
                Y_tilde(:,kk1) = reshape((Sel_mat*Y_tilde_k).',N_hv_sub*N_OFDM,1);
            end
            y_hat_UL = reshape(Y_tilde,K*N_hv_sub*N_OFDM,1);
            Y_bar_UL = reshape(y_hat_UL,N_OFDM,K*N_hv_sub).';

            N_miu_est_UL = TDU_ESPRIT_Algorithm(Y_bar_UL, N_U, K_miu_BS_Re, K_niu_BS_Re, K_miu_tau_Re, K_miu_BS_Im, K_niu_BS_Im, K_miu_tau_Im, M_sub, G_g, Sel_Mat_3D);
            if N_U == 1
                miu_BS_est = N_miu_est_UL(:,1);
                niu_BS_est = N_miu_est_UL(:,2);
                miu_tau_est = N_miu_est_UL(:,3);
            else
                [miu_tau_est, Index] = sort(N_miu_est_UL(:,3),'descend');
                niu_BS_est_temp = N_miu_est_UL(:,2);
                miu_BS_est_temp = N_miu_est_UL(:,1);
                niu_BS_est = niu_BS_est_temp(Index);
                miu_BS_est = miu_BS_est_temp(Index);
            end

            k_tau = (0:(K-1))';
            A_tau_est = exp(1i*k_tau*miu_tau_est.');
            m_BS = (0:N_bs(1)-1).';
            n_BS = (0:N_bs(2)-1).';
            A_miu_BS_est = exp(1i*m_BS*miu_BS_est.')/sqrt(N_bs(1));
            A_niu_BS_est = exp(1i*n_BS*niu_BS_est.')/sqrt(N_bs(2));
            A_BS_est = Khatri_Rao(A_niu_BS_est,A_miu_BS_est);
            A_temp = Khatri_Rao(Khatri_Rao(A_tau_est,(Sel_mat*W_tilde_UL'*A_BS_est)),S_pilot.');
            alpha_est = A_temp\y_hat_UL/sqrt(N_BS/N_U);

            H_f_est = zeros(N_BS,N_U,K);
            for kk_1 = 1:K
                H_f_est(:,:,kk_1) = A_BS_est*sqrt(N_BS/N_U)*diag(alpha_est.*exp(1i*(kk_1-1)*miu_tau_est));
            end

            difference_channel_MSE = zeros(K,1);
            true_channel_MSE = zeros(K,1);
            for kk_2 = 1:K
                difference_channel_MSE(kk_2) = norm(H_f_est(:,:,kk_2) - H_f(:,:,kk_2),'fro')^2;
                true_channel_MSE(kk_2) = norm(H_f(:,:,kk_2),'fro')^2;
            end
            NMSE_temp = sum(difference_channel_MSE)/sum(true_channel_MSE);
            NMSE(uu,ii) = NMSE(uu,ii) + NMSE_temp;

            if mod(iter,100) == 0
                toc
                disp(['  U = ' num2str(N_U_set(uu)) ', SNR = ' num2str(SNR_dBs(ii)) ', iter = ' num2str(iter)...
                    ' , NMSE_current = ' num2str(10*log10(NMSE_temp)) ' , NMSE_mean = ' num2str(10*log10(NMSE(uu,ii)/iter))])
            end
        end
        toc
        disp(['Finished ',num2str((uu-1)*length(SNR_dBs)+ii),'/', num2str(length(N_U_set)*length(SNR_dBs))]);
    end
end
disp('Finished all');
NMSE = NMSE/iterMax;
NMSE_dB_ESPRIT = 10*log10(NMSE);

%% Plot
MarkerSize = 6;
LineWidth = 2;
Fontsize = 15;

figure
plot(SNR_dBs,NMSE_dB_ESPRIT(1,:),'-.ko','MarkerSize',MarkerSize,'LineWidth',LineWidth);hold on;
plot(SNR_dBs,NMSE_dB_ESPRIT(2,:),'-.r^','MarkerSize',MarkerSize,'LineWidth',LineWidth);
plot(SNR_dBs,NMSE_dB_ESPRIT(3,:),'-.bs','MarkerSize',MarkerSize,'LineWidth',LineWidth);
axis([-10 15 -40 -5])
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('NMSE [dB]','Fontsize',Fontsize)
set(gca,'FontSize',Fontsize, 'linewidth',1.5);
set(gca, 'GridLineStyle', '-.'); grid on;
set(gcf, 'position', [700 300 650 550]);
h1 = legend('Proposed Scheme, {\it{U}} = 1', 'Proposed Scheme, {\it{U}} = 2', 'Proposed Scheme, {\it{U}} = 3','Location','southwest');
set(h1,'Fontsize',14);
