function [G0, G1] = linearize(fn_dyn, Xss, deriv)
% LINEARIZE find the coefficient matrices for a VAR(1).
%
%   Find the coefficient matrices for a linear system that can be written
%   as a VAR(1). Since the system is in the form, 
%       G0*X(t) = G1*X(t-1) + C + Psi*z(t) + Pi*eta(t),
%   we can find the G0 and G1 matrices by taking derivatives of the system
%   with respect to X(t) and X(t-1), holding the rest of the system at the
%   steady state, Xss.
%
%------------------------------------------------------------
%   INPUTS
%   - fn_dyn : function
%       Function governing model variable dynamics
%   - Xss : vector
%       Vector of model variables in steady state
%   - deriv : scalar
%       Scaling variable to determine step size of the derivative
%   OUTPUTS
%   - G0, G1 : matrices
%       First and second coefficient matrices in the linear system: 
%       G0*X(t) = G1*X(t-1) + C + Psi*z(t) + Pi*eta(t)
%------------------------------------------------------------


Fss = fn_dyn(Xss, Xss);

n_states = length(Xss);

% Find G0 matrix by taking derivatives of system with respect to X(t)
G0 = zeros(n_states,n_states);
for i = 1:n_states
	Xu = Xss;
	h = Xss(i) * deriv;
	if (h < deriv)
		h = deriv;
	end
	Xu(i) = Xu(i) + h;
	Fu = fn_dyn(Xu, Xss);
	G0(:,i) = (Fu - Fss)./h;
end

% Find G1 matrix by taking derivatives of system with respect to X(t-1)
G1 = zeros(n_states,n_states);
for i = 1:n_states
	Xu = Xss;
	h = Xss(i) * deriv;
	if (h < deriv)
		h = deriv;
	end
	Xu(i) = Xu(i) + h;
	Fu = fn_dyn(Xss, Xu);
	G1(:,i) = (Fu - Fss)./h;
end
G1 = -G1;


end