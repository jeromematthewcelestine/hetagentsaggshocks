function [Xss, idx] = investment_steadystate(opt, params)
%INVESTMENT_STEADYSTATE computes the steady state of the model
% 
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   OUTPUTS
%   - Xss : vector
%       Vector of model variables in steady state
%   - idx : structure
%       Contains the index values inside Xss that denote each
%       of the state variables' positions.
%       
%---------------------------------

kp_grid = repmat(opt.k_grid, 1, opt.n_z);

x_ss = 1;
Ex_ss = x_ss;
price_ss = 1;
Eprice_ss = price_ss;

converged = 0;
for iter = 1:opt.n_iter
	kp_grid_old = kp_grid;
	
	kp_grid = investment_solve_for_policy(opt, params, kp_grid_old, Ex_ss, price_ss, Eprice_ss);

	error = max(abs(kp_grid(:) - kp_grid_old(:)));
	if (error < opt.kp_error_tol)
		converged = 1;
		break;
	end
end
kp_ss = kp_grid;

fprintf('converged %d, iter %d, error %f\n',converged, iter, error);

Q_ss	= compute_transition_matrix(opt,kp_ss);
dist_ss = compute_stationary_distribution(opt, Q_ss);

aggs_ss = compute_aggregates(opt, params, dist_ss, kp_ss, x_ss);

idx.kp			= 1:opt.n_kp;
idx.Ekp			= opt.n_kp+1:2*opt.n_kp;
idx.dist		= 2*opt.n_kp+1:3*opt.n_kp;
idx.x			= 3*opt.n_kp+1;
idx.Ex			= 3*opt.n_kp+2;
idx.output		= 3*opt.n_kp+3;
idx.investment	= 3*opt.n_kp+4;
idx.consumption	= 3*opt.n_kp+5;
idx.price		= 3*opt.n_kp+6;
idx.Eprice		= 3*opt.n_kp+7;

n_states = idx.Eprice;

Xss = zeros(n_states,1);
Xss(idx.kp)				= reshape(kp_ss,opt.n_kp,1);
Xss(idx.Ekp)			= reshape(kp_ss,opt.n_kp,1);
Xss(idx.dist)			= dist_ss;
Xss(idx.x)				= log(x_ss);
Xss(idx.Ex)				= log(x_ss);
Xss(idx.price)			= log(price_ss);
Xss(idx.Eprice)			= log(price_ss);
Xss(idx.consumption)	= log(aggs_ss.consumption);
Xss(idx.output)			= log(aggs_ss.output);
Xss(idx.investment)		= log(aggs_ss.investment);


end
