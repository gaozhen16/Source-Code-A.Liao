clear,clc

N_BS = 64;
N_MS = 64;
n_BS = (0:N_BS-1).';
n_MS = (0:N_MS-1).';
N_RF_BS = 4;
N_RF_MS = 4;
fc = 30e9;
lambda = 3e8/fc;
d_ant = lambda/2;
Lp = 5;
Ns = 3;
m1 = 13; m2 = 12;
N = 10;

[F_Mat, F, W_H_Mat, W_H, Lambda] = Hybrid_codebook_Com(N_BS, N_MS, N_RF_BS, N_RF_MS, N, Ns);
W_H_blkdiag = zeros(N*Ns,N*N_BS);
for n1 = 1:N
    W_H_blkdiag((n1-1)*Ns+1:n1*Ns,(n1-1)*N_BS+1:n1*N_BS) = W_H_Mat(:,:,n1);
end
S_pilot = dftmtx(Ns);
S_blkdiag = kron(eye(N),S_pilot);

iterMax = 1000;
PNR_dBs = -5:5:20;
NMSE_ESPRIT_L5 = zeros(1,length(PNR_dBs));

for i_PNR_dB = 1:length(PNR_dBs)
    tic
    sigma2 = 10^(-(PNR_dBs(i_PNR_dB)/10));
    sigma = sqrt(sigma2);
    for iter = 1:iterMax
        [H_up, A_BS, A_MS, D, theta, phi] = mmWave_uplink_channel(N_BS, N_MS, Lp, lambda, d_ant);
        Y_Mat_com = W_H*H_up*F*S_blkdiag + sigma*W_H_blkdiag*(randn(N*N_BS,N*Ns) + 1i*randn(N*N_BS,N*Ns))/sqrt(2);
        H_NN_hat_com = Y_Mat_com*S_blkdiag'/Ns;
        [theta_est, phi_est] = TD_esprit_theta_phi(H_NN_hat_com, m1, m2, Lp, lambda, d_ant);
        A_BS_ESPRIT = exp(1i*2*pi/lambda*d_ant*n_BS*sin(theta_est).')/sqrt(N_BS);
        A_MS_ESPRIT = exp(1i*2*pi/lambda*d_ant*n_MS*sin(phi_est).')/sqrt(N_MS);
        [H_hatm, H_hatn] = size(H_NN_hat_com); H_NN_hat_com_vec = reshape(H_NN_hat_com,H_hatm*H_hatn,1);
        A_mat_l = W_H*A_BS_ESPRIT; A_mat_r = A_MS_ESPRIT'*F;
        A_mat = Khatri_Rao((A_mat_r).', A_mat_l);
        D_ESPRIT_col = A_mat\H_NN_hat_com_vec;
        H_up_ESRIT = A_BS_ESPRIT*diag(D_ESPRIT_col)*A_MS_ESPRIT';
        NMSE_ESPRIT_temp = norm(H_up_ESRIT - H_up,'fro')^2/norm(H_up,'fro')^2;
        NMSE_ESPRIT_L5(i_PNR_dB) = NMSE_ESPRIT_L5(i_PNR_dB) + NMSE_ESPRIT_temp;
        
        if mod(iter,200) == 0
            disp(['  snr = ' num2str(PNR_dBs(i_PNR_dB)) ', iter = ' num2str(iter)  ' , NMSE_ESPRIT = ' num2str(NMSE_ESPRIT_L5(i_PNR_dB)/iter)])
        end
    end
    NMSE_ESPRIT_L5(i_PNR_dB) = NMSE_ESPRIT_L5(i_PNR_dB)/iterMax;
    toc
    disp(['Finished ',num2str(i_PNR_dB),'/', num2str(length(PNR_dBs)) ' , NMSE_ESPRIT = ' num2str(NMSE_ESPRIT_L5(i_PNR_dB))]);
end
NMSE_ESPRIT_dB_L5 = 10*log10(NMSE_ESPRIT_L5);

disp('Finished all');


%% Plot
figure
plot(PNR_dBs,NMSE_ESPRIT_dB_L5,'-ko','LineWidth',1.5); grid on;
xlabel('SNR [dB]'),ylabel('NMSE [dB]')
