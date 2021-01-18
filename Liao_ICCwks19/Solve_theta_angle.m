function theta_angle = Solve_theta_angle(Phi_com, N_miu, p, q)

c_vec_com = zeros(5,N_miu);
for n_miu = 1:N_miu
	Phi_n = Phi_com(:,:,n_miu);
	if q - p > 1
        c_vec_pq_n = c_vec_pq(Phi_n, p, q);
        c_add_temp = zeros(5,q-p-1);
        for kk = p+1:q-1
            v_pk_n = Phi_n(p,kk);
            v_qk_n = Phi_n(q,kk);
            v_kp_n = Phi_n(kk,p);
            v_kq_n = Phi_n(kk,q);
            c_add_temp(:,kk-p) = c_vec_add(v_pk_n,v_qk_n) - c_vec_add(v_kp_n,v_kq_n);
        end
        c_add_temp_sum = sum(c_add_temp,2);
        c_vec_n = c_vec_pq_n + c_add_temp_sum;
	else
        c_vec_n = c_vec_pq(Phi_n, p, q);
	end
	c_vec_com(:,n_miu) = c_vec_n;
end
c_vec_sum = sum(c_vec_com,2);
pt_roots = roots((flipud(eye(5))*c_vec_sum).');
real_roots_pt = real(pt_roots(abs(imag(pt_roots))<1e-18));
if isempty(real_roots_pt)
    theta_angle = 0;
else
	num_real_roots = length(real_roots_pt);
	derivative_pt = zeros(num_real_roots,1);
	for i_num = 1:num_real_roots
        t = real_roots_pt(i_num);
        derivative_pt(i_num) = c_vec_sum(2)+2*c_vec_sum(3)*t+3*c_vec_sum(4)*t^2+4*c_vec_sum(5)*t^3;
	end
	real_roots_pt_used = real_roots_pt(derivative_pt>0);
	if isempty(real_roots_pt_used)
        theta_angle = 0;
    else
        theta_angle = atan(real_roots_pt_used);
	end
end

end

function c_vec_pq_n = c_vec_pq(Phi_n, p, q)
c_vec_pq_n = zeros(5,1);
v_pp_n = Phi_n(p,p);
v_pq_n = Phi_n(p,q);
v_qp_n = Phi_n(q,p);
v_qq_n = Phi_n(q,q);
c_vec_pq_n(1) =v_qp_n*(v_pp_n-v_qq_n);
c_vec_pq_n(2) = (v_pp_n-v_qq_n)^2-2*v_qp_n*(v_pq_n+v_qp_n);
c_vec_pq_n(3) = -3*(v_pp_n-v_qq_n)*(v_pq_n+v_qp_n);
c_vec_pq_n(4) = -(v_pp_n-v_qq_n)^2+2*v_pq_n*(v_pq_n+v_qp_n);
c_vec_pq_n(5) = v_pq_n*(v_pp_n-v_qq_n);
end

function c_vec_add = c_vec_add(a, b)
c_vec_add = zeros(5,1);
c_vec_add(1) = a*b;
c_vec_add(2) = a^2-b^2;
c_vec_add(3) = 0;
c_vec_add(4) = a^2-b^2;
c_vec_add(5) = -a*b;
end