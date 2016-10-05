function dist = make_stationary_distribution(opt, transition_matrix)

n_dist = opt.n_k * opt.n_z;
dist = 1./n_dist * ones(n_dist,1);

for dist_iter = 1:opt.n_dist_iter
	dist_old	= dist;
	dist		= transition_matrix' * dist_old;
	
	error = max(abs(dist(:) - dist_old(:)));
	if (error < opt.dist_error_tol)
		converged = 1;
		break;
	end
end
