function dist = compute_stationary_distribution(opt, transition_matrix)
%COMPUTE_STATIONARY_DISTRIBUTION stationary distribution of the model
%
%   Compute the stationary distribution of the model. This does not feature
%   aggregate uncertainty, although it does contain idiosyncratic
%   uncertainty. 
%   
%   Compute by iterating on the discretized density transition equation:
%       dist(t) = Q(t)*dist(t-1)
%   where Q(t) = transition_matrix. 
%
%------------------------------------------------------------
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - transition_matrix : matrix
%       Transition probability matrix from points in the state space to 
%       points in the state space next period.
%   
%   OUTPUTS
%   - dist : vector
%       
%------------------------------------------------------------

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


end