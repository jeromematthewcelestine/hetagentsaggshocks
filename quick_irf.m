function state = quick_irf(A, B, T)
% QUICK_IRF plots impulse response functions for the model
%
%   INPUTS
%   - A : matrix
%     VAR(1) coefficient matrix  
%   - B : matrix
%     Impact coefficient matrix on unit shocks to the VAR(1) system
%   - T : scalar
%     Length of the IRF.
%   OUTPUTS
%   - state : matrix
%     Impulse responses 
%------------------------------------------------------------

state = zeros(size(A,1), T);
innov = zeros(size(B,2), T);
innov(1,2) = 0.01;

for t = 2:T
	state(:,t) = A * state(:,t-1) + B * innov(:,t);
end



end
