function Q_N = Q(N)
% Generating a unitary matrix Q_N, which defined as follows:
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