function [G0, G1] = linearize(fn_dyn, Xss, deriv)
% LINEARIZE finds the coefficient matrices in the linear VAR(1).
%
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
%       G0*x(t) = G1*x(t-1) + C + Psi*z(t) + Pi*eta(t)
%
%------------------------------------------------------------


Fss = fn_dyn(Xss, Xss);

n_states = length(Xss);

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