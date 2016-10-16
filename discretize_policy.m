function Q = discretize_policy(grid, pol)
%DISCRETIZE_POLICY compute transition probability matrix
%
%   Creates a transition probability matrix using linear interpolation over
%   the policy function.
%
%------------------------------------------------------------
%   INPUTS
%   - grid : vector
%       Define grid of points for the given state variable
%   - pol : vector
%       Policy function.
%   OUTPUTS
%   - Q : matrix
%       Transition probability matrix from policy function to points on the
%       variable's grid (i.e. from k' to kgrid).
%------------------------------------------------------------

n_grid = length(grid);

Q = zeros(length(pol),length(grid));
for node_idx = 1:length(pol)
	left_idx = find(pol(node_idx) > grid, 1, 'last');
	if (isempty(left_idx))
		Q(node_idx,1) = 1;
	elseif (left_idx == n_grid)
		Q(node_idx,n_grid) = 1;
	else
		Q(node_idx,left_idx) = (grid(left_idx+1) - pol(node_idx))/(grid(left_idx+1) - grid(left_idx));
		Q(node_idx,left_idx+1) = 1 - Q(node_idx,left_idx);
	end
end



end