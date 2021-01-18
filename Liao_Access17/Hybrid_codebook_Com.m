function [F_Mat, F, W_H_Mat, W_H, Lambda] = Hybrid_codebook_Com(N_BS, N_MS, N_RF_BS, N_RF_MS, N, Ns)

    if N == 1
        error('N must be greater than 1');
    end
    if N*Ns <= N_MS && N*Ns <= N_BS
        N_tr = N;
    else
        error('N is too big');
    end

    F_Mat = zeros(N_MS,Ns,N_tr);
    F_dft_MS = dftmtx(N_RF_MS)/sqrt(N_MS);
    F_BB = F_dft_MS(:,1:Ns)*sqrt(N_MS);
    for n1 = 1:N_tr
        F_RF_temp = repmat(F_dft_MS(:,N_RF_MS)',N_MS,1);
        F_RF_temp((n1-1)*Ns+1:n1*Ns,:) = F_dft_MS(:,1:Ns)';
        F_Mat(:,:,n1) = F_RF_temp*F_BB;
    end
    Lambda = sqrt(N_RF_MS)/norm(F_Mat(:,:,1),'fro');    % Ns
    F_Mat = Lambda*F_Mat;
    F = [];
    for n2 = 1:N_tr
        F_temp = F_Mat(:,:,n2);
        F = [F,F_temp];
    end

    W_H_Mat = zeros(Ns,N_BS,N_tr);
    F_dft_BS = dftmtx(N_RF_BS)/sqrt(N_BS);
    W_BB = F_dft_BS(:,1:Ns);
    for n3 = 1:N_tr
        W_RF_temp = repmat(F_dft_BS(:,N_RF_BS)',N_BS,1);
        W_RF_temp((n3-1)*Ns+1:n3*Ns,:) = F_dft_BS(:,1:Ns)';
        W_H_Mat(:,:,n3) = (W_RF_temp*W_BB)';
    end
    W_H = [];
    for n4 = 1:N_tr
        W_H_temp = W_H_Mat(:,:,n4);
        W_H = [W_H;W_H_temp];
    end

end