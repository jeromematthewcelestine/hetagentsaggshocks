function [G0, G1] = linearize(fn_dyn, Xss, deriv)

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