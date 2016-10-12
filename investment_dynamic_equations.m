function f = investment_dynamic_equations(opt, params, X_t, X_tm1, Xss, idx)
%INVESTMENT_DYNAMIC_EQUATIONS compute change in variables
%
%   Compute the difference between the LHS and RHS of a variable. That is,
%   compare the given value of a variable at t with the value of the
%   variable if constructed from all other variables in the system. The
%   function can then be used to compute derivatives by varying various
%   variables in the system.
% 
%---------------------------------
%   INPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%   - X_t : vector
%       The vector of state variables in the current period
%   - X_tm1 : vector
%       The vector of state variables in previous period (i.e. "X at t
%       minus 1")
%   - Xss : vector
%       Vector of model variables in steady state
%   - idx : structure
%       Contains the index values inside Xss that denote each
%       of the state variables' positions. 
% 
%   OUTPUTS
%   - f : vector
%       For every variable in the model, compute the difference between the
%       variable and the construction of that variable using the other
%       variables. E.g. suppose x_t = y_t + z_t, then f = x_t - (y_t + z_t)
%---------------------------------

n_states	= length(Xss);

% Extract model variables
kp_grid_t			= X_t(idx.kp);
Ekp_grid_t			= X_t(idx.Ekp);
dist_t				= X_t(idx.dist);
x_t					= X_t(idx.x);
Ex_t				= X_t(idx.Ex);
output_t			= exp(X_t(idx.output));
investment_t		= exp(X_t(idx.investment));
consumption_t		= exp(X_t(idx.consumption));
price_t				= exp(X_t(idx.price));
Eprice_t			= exp(X_t(idx.Eprice));

kp_grid_tm1			= X_tm1(idx.kp);
Ekp_grid_tm1		= X_tm1(idx.Ekp);
dist_tm1			= X_tm1(idx.dist);
x_tm1				= X_tm1(idx.x);
Ex_tm1				= X_tm1(idx.Ex);
output_tm1			= exp(X_tm1(idx.output));
investment_tm1		= exp(X_tm1(idx.investment));
consumption_tm1		= exp(X_tm1(idx.consumption));
price_tm1			= exp(X_tm1(idx.price));
Eprice_tm1			= exp(X_tm1(idx.Eprice));


% Compute capital policy functions
kp_grid_new			= investment_solve_for_policy(opt, params, reshape(Ekp_grid_t,opt.n_k,opt.n_z), exp(Ex_t), price_t, Eprice_t);
kp_grid_new			= reshape(kp_grid_new, opt.n_kp, 1);

% Compute density
% NOTE: dist_tm1 is the distribution over states in time t that occured as
% a result of decisions made at time t-1 (i.e. k_t is the result of capital
% choices made at t-1). 
new_Q				= compute_transition_matrix(opt, kp_grid_t);
dist_new			= new_Q' * dist_tm1;

% Compute aggregate output, investment, and consumption
output_mesh			= exp(x_t) * opt.z_mesh .* (opt.k_mesh.^params.alpha);
investment_mesh		= reshape(kp_grid_new,opt.n_k,opt.n_z) - (1-params.delta).*opt.k_mesh;

output_new			= sum(sum( output_mesh .* reshape(dist_tm1,opt.n_k,opt.n_z) ));
investment_new		= sum(sum( investment_mesh .* reshape(dist_tm1,opt.n_k,opt.n_z) ));
% consumption_new		= output_new - investment_new;
consumption_new		= output_t - investment_t;

% Compute prices
price_new			= 1/consumption_t;

% Compute aggregate state variable
x_new				= params.rho * x_tm1;

% Compute differences, f, and pack up output
f = zeros(n_states,1);
f(idx.kp)			= kp_grid_t				- kp_grid_new;
f(idx.Ekp)			= kp_grid_t				- Ekp_grid_tm1;
f(idx.dist)			= dist_t				- dist_new;
f(idx.x)			= x_t					- x_new;
f(idx.Ex)			= x_t					- Ex_tm1;
f(idx.output)		= log(output_t)			- log(output_new);
f(idx.investment)	= log(investment_t)		- log(investment_new);
f(idx.consumption)	= log(consumption_t)	- log(consumption_new);
f(idx.price)		= log(price_t)			- log(price_new);
f(idx.Eprice)		= log(price_t)			- log(Eprice_tm1);

end

