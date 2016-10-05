clear;
clc;

tic;

opt.n_iter			= 100;
opt.n_dist_iter		= 1000;
opt.kp_error_tol	= 1e-10;
opt.dist_error_tol	= 1e-8;
opt.derivative_size = 1e-5;

params.beta		= 0.9;
params.alpha	= 0.2;
params.delta	= 0.1;
params.rho		= 0.9;
params.phi		= 0.1;

optimality_tol = 1e-15;
opt.fsolve_options = optimoptions('fsolve','Display','none','OptimalityTolerance',optimality_tol);

error_tol = 1e-10;
opt.k_min = 0.01;
opt.k_max = 1;
opt.n_k = 10;
opt.k_grid			= exp(linspace(log(opt.k_min),log(opt.k_max),opt.n_k))';
opt.k_grid(1)		= opt.k_min;
opt.k_grid(opt.n_k) = opt.k_max;

opt.n_z = 2;
opt.z_grid = [0.7, 0.9];
opt.Pz = [0.85, 0.15; 0.05, 0.95];
opt.Qz = kron(opt.Pz, ones(opt.n_k,1));

opt.n_kp = opt.n_z * opt.n_k;

opt.z_mesh	= repmat(opt.z_grid,opt.n_k,1);
opt.k_mesh	= repmat(opt.k_grid,1,opt.n_z);

addpath('Het1');

[Xss, idx] = Investment_steadystate(opt, params);

%%

f_bnd = @(X_t, X_tm1) Investment_dynamic_equations(opt, params, X_t, X_tm1, Xss, idx);

%% 
n_states = length(Xss);

Psi_r1 = zeros(n_states,1);
Psi_r1(idx.x) = 1;

Pi_r1 = zeros(n_states,opt.n_kp+2);
Pi_r1(idx.Ekp,1:opt.n_kp) = eye(opt.n_kp);
Pi_r1(idx.Ex,opt.n_kp+1) = 1;
Pi_r1(idx.Eprice,opt.n_kp+2) = 1;

C_r1 = zeros(n_states,1);

[G0_r1, G1_r1] = linearize(f_bnd, Xss, opt.derivative_size);

%% 

[A_r1, ~, B_r1, ~, ~, ~, ~, eu] = gensys(G0_r1, G1_r1, C_r1, Psi_r1, Pi_r1);
eu

toc;

%% plotting

% plot ss
kp_grid_ss	= reshape(Xss(idx.kp),opt.n_k,opt.n_z);
hist_ss		= reshape(Xss(idx.dist),opt.n_k,opt.n_z);

figure;
plot(opt.k_grid,kp_grid_ss)
figure;
plot(opt.k_grid,hist_ss)

% plot IRFs
figure;
state = quick_irf(A_r1, B_r1);
plot(state(idx.x:idx.Eprice,:)')
legend('x','Ex','y','i','c','p','Ep')

