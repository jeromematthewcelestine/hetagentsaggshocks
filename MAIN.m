%---------------------------------
% 
% Demonstration of the Reiter method for solving a model with heterogeneous 
% agents and aggregate shocks in general equilibrium. Solves a simple model 
% of firm investment with persistent aggregate and idiosyncratic productivity 
% shocks.
%
% Note, this code uses the Gensys code of Chris Sims. This can be
% downloaded from: 
%   http://sims.princeton.edu/yftp/gensys/
% 
% Reiter (2009), Journal of Economic Dynamics & Control, 
% "Solving heterogeneous-agent models by projection and perturbation."
%
% 
% Written by Jerome Williams (2016)
% Edited by James Graham.
%---------------------------------


%% Setup 
close all;
clear;
clc;

tic;

[opt, params] = setup;


%% Compute steady state

[Xss, idx] = investment_steadystate(opt, params);


%% Construct linear system matrices
% Input the system in the form: G0_r1*X(t) = G1_r1*X(t-1) + C_r1 + Psi_r1*z(t) + Pi_r1*eta(t)

n_states = length(Xss);

Psi_r1 = zeros(n_states,1);
Psi_r1(idx.x) = 1;

Pi_r1 = zeros(n_states,opt.n_kp+2);
Pi_r1(idx.Ekp,1:opt.n_kp) = eye(opt.n_kp);
Pi_r1(idx.Ex,opt.n_kp+1) = 1;
Pi_r1(idx.Eprice,opt.n_kp+2) = 1;

C_r1 = zeros(n_states,1);

f_bnd = @(X_t, X_tm1) investment_dynamic_equations(opt, params, X_t, X_tm1, Xss, idx);
[G0_r1, G1_r1] = linearize(f_bnd, Xss, opt.derivative_size);


%% Use Gensys to solve the linear Rational Expectations system 
%
% Returns system of the form: X(t) = A_r1*X(t-1) + Ctilde + B_r1*z(t) + ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% Note, here Ctilde=0, and we are interested in A_r1 and B_r1 

[A_r1, ~, B_r1, ~, ~, ~, ~, eu] = gensys(G0_r1, G1_r1, C_r1, Psi_r1, Pi_r1);

fprintf('--------------------------------\n')
fprintf('Existence of eqm: \t %i\n', eu(1))
fprintf('Uniqueness of eqm: \t %i\n', eu(2))
fprintf('--------------------------------\n')

toc;

%% Plot 

fontsize = 12;

% Plot steady state policy functions and distribution
kp_grid_ss	= reshape(Xss(idx.kp),opt.n_k,opt.n_z);
hist_ss		= reshape(Xss(idx.dist),opt.n_k,opt.n_z);

figure;
subplot(1,2,1)
plot(opt.k_grid,kp_grid_ss, 'linewidth', 2)
grid on
title('Policy function')
xlabel('Current capital, k')
legend('Low productivity','High productivity')
set(gca, 'fontsize', fontsize)

subplot(1,2,2)
plot(opt.k_grid,hist_ss, 'linewidth', 2)
grid on
title('Distribution of capital in SS')
xlabel('Current capital, k')
legend('Low productivity','High productivity')
set(gca, 'fontsize', fontsize)


% Plot impulse response functions to aggregate shock
T = 50;
state = quick_irf(A_r1, B_r1, T);

figure;
plot(state(idx.x:idx.Eprice,:)')
grid on
title('Impulse response functions')
legend('x','Ex','y','i','c','p','Ep')
set(gca, 'fontsize', fontsize)
