function kp_grid_out = investment_solve_for_policy(opt, params, kp_grid_in, Ex, price, Eprice)
%INVESTMENT_SOLVE_FOR_POLICY solve for individual firm policy functions
% 
%   Given the firm's objective function and a current guess for the policy
%   function, update the guess for the firm's optimal capital choice. In
%   this case, the firm's objective function is defined as the difference
%   between the left and right hand side of the FOC. Hence, we are solving
%   for the zero of that system.
%
%---------------------------------
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   - kp_grid_in : vector
%       Current guess for the capital policy.
%   - Ex : scalar
%       Expected value of next period's aggregate state.
%   - price : scalar
%       Current price of the investment good.
%   - Eprice : scalar
%       Expected value of next period's price of the investment good.
% 
%   OUTPUTS
%   - kp_grid_out : vector
%       Updated guess for the capital policy.
%       
%---------------------------------

kp_grid_out = zeros(size(kp_grid_in));

for k_idx = 1:opt.n_k
	for z_idx = 1:opt.n_z
		k = opt.k_grid(k_idx);
		
		obj_fn = @(kp) investment_firm_objective_function(opt, params, kp_grid_in, ...
                                                          k, z_idx, kp, Ex, price, Eprice);

		% 	first check if lower bound binds
		diff_k = obj_fn(opt.k_grid(1));

		if (diff_k > 0)
			% lower bound binds
			kp_star		= opt.k_min;
		else
			kp_star = fsolve(obj_fn, kp_grid_in(k_idx,z_idx), opt.fsolve_options);
		end
		kp_grid_out(k_idx,z_idx)	= kp_star;
	end
end


end