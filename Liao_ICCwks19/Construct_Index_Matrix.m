function Index_Mat = Construct_Index_Matrix(Nr, N_sub, Ns)

N_H = Nr(1);
N_V = Nr(2);
N_h_sub = N_sub(1);
N_v_sub = N_sub(2);
Block_Matrix = [eye(N_h_sub);zeros(N_H-N_h_sub,N_h_sub)];
N_b = ceil(N_h_sub*N_v_sub/Ns);
Sel_Matrix = [eye(Ns*N_b);zeros((N_v_sub+1)*N_h_sub-Ns*N_b,Ns*N_b)];
Index_Mat = [kron(eye(N_v_sub+1),Block_Matrix);zeros(N_V*N_H-(N_v_sub+1)*N_H,(N_v_sub+1)*N_h_sub)]*Sel_Matrix;

end