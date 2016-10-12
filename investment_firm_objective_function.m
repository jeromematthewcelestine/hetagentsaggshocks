function ret = investment_firm_objective_function(opt, params, kpp_grid, k, z_idx, kp, x_tp1, price_t, price_tp1)
%INVESTMENT_FIRM_OBJECTIVE_FUNCTION 
% 
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   - kpp_grid
%   
%   - k
% 
%   - z_idx
%
%   - kp
% 
%   - x_tp1
%
%   - price_t 
%
%   - price_tp1
%
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
%   - ret : 
%       
%       
%---------------------------------


lhs = price_t * (1 + 2*params.phi*(kp./k) + 2*params.phi*(1-params.delta));

kpp = [interp1(opt.k_grid,kpp_grid(:,1),kp), interp1(opt.k_grid,kpp_grid(:,2),kp)];
rhs = price_tp1 * params.beta * opt.Pz(z_idx,:) * ...
	(x_tp1 .* opt.z_grid' * params.alpha * kp.^(params.alpha-1) + 1 - params.delta - ...
                          params.phi * (kpp'./kp).^2 + params.phi * (1-params.delta)^2);
ret = lhs - rhs;

end