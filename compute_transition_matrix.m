function Q = compute_transition_matrix(opt, kp_grid)
%COMPUTE_TRANSITION_MATRIX creates transition probaility matrix
%
%
%------------------------------------------------------------
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - kp_grid : vector
%       Define grid of points for the capital state
%   OUTPUTS
%   - Q : matrix
%       Transition probability matrix from [pints in the state space to 
%       points in the state space next period.
%------------------------------------------------------------


Qk = discretize_policy(opt.k_grid, reshape(kp_grid, opt.n_k*opt.n_z, 1));
Q = repelem(opt.Qz,1,opt.n_k).*repmat(Qk,1,opt.n_z);

end