function state = quick_irf(A, B)

T = 50;
state = zeros(size(A,1), T);
innov = zeros(size(B,2), T);
innov(1,2) = 0.01;
for t = 2:T
	state(:,t) = A * state(:,t-1) + B * innov(:,t);
end
