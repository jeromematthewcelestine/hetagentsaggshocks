function [ opt, params ] = setup()
%SETUP Runs the setup for the Reiter method. 
% 
%   Setup the options and parameters for the model. 
%---------------------------------
%   INPUTS
%   - 
%   OUTPUTS
%   - opt : structure
%       Options for number of iterations and tolerances on approximation
%       routines
%   - params : structure
%       Model parameters
%
%---------------------------------

%% Parameters

params.beta         = 0.9;
params.alpha        = 0.2;
params.delta        = 0.1;
params.rho          = 0.9;
params.phi          = 0.1;


%% Options
opt.n_iter			= 100;
opt.n_dist_iter		= 1000;
opt.kp_error_tol	= 1e-10;
opt.dist_error_tol	= 1e-8;
opt.derivative_size = 1e-5;

optimality_tol      = 1e-15;
opt.fsolve_options  = optimoptions('fsolve','Display','none','OptimalityTolerance',optimality_tol);

error_tol           = 1e-10;
opt.k_min           = 0.01;
opt.k_max           = 1;
opt.n_k             = 25;
opt.k_grid			= exp(linspace(log(opt.k_min),log(opt.k_max),opt.n_k))';
opt.k_grid(1)		= opt.k_min;
opt.k_grid(opt.n_k) = opt.k_max;

% Idiosyncratic state grids and Markov transition matrices.
opt.n_z             = 2;            
opt.z_grid          = [0.7, 0.9];
opt.Pz              = [0.85, 0.15; 0.05, 0.95];
opt.Qz              = kron(opt.Pz, ones(opt.n_k,1));

opt.n_kp            = opt.n_z * opt.n_k;

opt.z_mesh          = repmat(opt.z_grid,opt.n_k,1);
opt.k_mesh          = repmat(opt.k_grid,1,opt.n_z);


end

