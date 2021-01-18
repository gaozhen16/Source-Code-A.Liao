function [W_tilde, W_Mat] = Combiner_Design(Nr, N_sub, N_RF, Ns, N_bits)

    N_R = Nr(1)*Nr(2);
    N_h_sub = N_sub(1);
    N_v_sub = N_sub(2);
    if N_h_sub > Nr(1) || N_v_sub > Nr(2)
        error('M_h_sub or M_v_sub is wrong!');
    end
    Index_Matrix = Construct_Index_Matrix(Nr, N_sub, Ns);
    if N_RF == 2 || log2(N_RF) == fix(log2(N_RF))
        Orth_mtx = hadamard(N_RF);
    else
        Orth_mtx = dftmtx(N_RF);
    end
    W_BB_k = Orth_mtx(:,1:Ns);
    N_b = ceil(N_h_sub*N_v_sub/Ns);
    W_Mat = zeros(N_R,Ns,N_b);
    W_tilde = [];
    for n_b = 1:N_b
        Index_Sub_Matrix = Index_Matrix(:,(n_b-1)*Ns+1:n_b*Ns);
        index = mod(find(Index_Sub_Matrix~=0),N_R);
        W_RF_temp = repmat(Orth_mtx(:,N_RF)',N_R,1);
        W_RF_temp(index.',:) = W_BB_k';
        W_RF = Quantize(W_RF_temp, N_bits);
        W = real(W_RF*W_BB_k);
        W_Mat(:,:,n_b) = W;
        W_tilde = [W_tilde, W];
    end
end

function Quantized_Matrix = Quantize(Matix, N_bits)
% Quantizing the analog precoding/conmbining matrix

N_Bits = 2^N_bits;
N_R = size(Matix,1);
Matix_quan_phase = floor((angle(Matix)+pi)*N_Bits/(2*pi)) *2*pi/N_Bits - pi;
Quantized_Matrix = exp(1i*Matix_quan_phase)/sqrt(N_R);

end