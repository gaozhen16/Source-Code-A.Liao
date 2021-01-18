function Index_Matrix = Construct_Index_Matrix(Nr_or_Nt, N_sub, Ns)

    Mrt = Nr_or_Nt(1);
    Nrt = Nr_or_Nt(2);
    Mrt_sub = N_sub(1);
    Nrt_sub = N_sub(2);
    Block_Matrix = [eye(Mrt_sub);zeros(Mrt-Mrt_sub,Mrt_sub)];
    N_b = ceil(Mrt_sub*Nrt_sub/Ns);
    Sel_Matrix = [eye(Ns*N_b);zeros((Nrt_sub+1)*Mrt_sub-Ns*N_b,Ns*N_b)];
    Index_Matrix = [kron(eye(Nrt_sub+1),Block_Matrix);zeros(Nrt*Mrt-(Nrt_sub+1)*Mrt,(Nrt_sub+1)*Mrt_sub)]*Sel_Matrix;

end