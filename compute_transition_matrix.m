function Q = compute_transition_matrix(opt, kp_grid)

Qk = discretize_policy(opt.k_grid, reshape(kp_grid, opt.n_k*opt.n_z, 1));
Q = repelem(opt.Qz,1,opt.n_k).*repmat(Qk,1,opt.n_z);