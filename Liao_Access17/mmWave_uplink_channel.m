function [H_uplink, A_BS, A_MS, D, theta, phi] = mmWave_uplink_channel(N_BS, N_MS, Lp, lambda, d_ant)

ang_min = -pi/3;ang_max = pi/3;
n_BS = (0:N_BS-1).';
n_MS = (0:N_MS-1).';
alpha = sort((randn(Lp,1)+1i*randn(Lp,1))/sqrt(2),'descend');
D = diag(sqrt(N_BS*N_MS/Lp)*alpha);
theta = ang_min+(ang_max-ang_min).*rand(Lp,1);
phi = ang_min+(ang_max-ang_min).*rand(Lp,1);
A_BS = exp(1i*2*pi/lambda*d_ant*n_BS*sin(theta.'))/sqrt(N_BS);
A_MS = exp(1i*2*pi/lambda*d_ant*n_MS*sin(phi.'))/sqrt(N_MS);

H_uplink = A_BS*D*A_MS';

end