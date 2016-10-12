function kp_grid_out = investment_solve_for_policy(opt, params, kp_grid_in, Ex, price, Eprice)
%INVESTMENT_SOLVE_FOR_POLICY solve for individual firm policy functions
% 
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   - kp_grid_in : 
% 
%   - Ex : 
%       
%   - price : 
%      
%   - Eprice : 
% 
% 
%   OUTPUTS
%   - kp_grid_out : vector
%       
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