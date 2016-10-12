function ret = investment_firm_objective_function(opt, params, kpp_grid, k, z_idx, kp, x_tp1, price_t, price_tp1)
%INVESTMENT_FIRM_OBJECTIVE_FUNCTION compute value of the objective function
% 
%   Compute value of the objective function given current and future
%   prices, future aggregate state, current idiosyncratic state variables,
%   and the current and next period's capital policy.
%
%   In this case, the firm's objective function is defined as the difference
%   between the left and right hands side of the first order condition. At
%   an optimum, this difference is zero. 
%---------------------------------
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   - kpp_grid : vector
%       Current guess for the capital policy function next period.
%   - k : scalar
%       Current value of the capital state.
%   - z_idx : scalar
%       Index number of the idiosyncratic productivity state.
%   - kp : scalar
%       Current guess for the capital policy in the current (k,z) state.
%   - x_tp1 : scalar
%       Value of next period's aggregate state.
%   - price_t : scalar
%       Current price of the investment good.
%   - price_tp1 : scalar
%       Value of next period's price of the investment good.
% 
%   OUTPUTS
%   - ret : 
%       Difference between the left and right hand side of the FOC.
%       
%---------------------------------


lhs = price_t * (1 + 2*params.phi*(kp./k) + 2*params.phi*(1-params.delta));

kpp = [interp1(opt.k_grid,kpp_grid(:,1),kp), interp1(opt.k_grid,kpp_grid(:,2),kp)];

rhs = price_tp1 * params.beta * opt.Pz(z_idx,:) * ...
	(x_tp1 .* opt.z_grid' * params.alpha * kp.^(params.alpha-1) + 1 - params.delta - ...
                          params.phi * (kpp'./kp).^2 + params.phi * (1-params.delta)^2);

ret = lhs - rhs;

end