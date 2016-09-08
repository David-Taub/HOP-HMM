
% v - 1 x k - vector of indices
% M - n x k
% example:
% [3,3,2,1], 3 -> 
% 0 0 0 1
% 0 0 1 0
% 1 1 0 0 
function M = vec2mat(v, n)
	M = full(sparse(v, 1:length(v), 1, n, length(v))) == 1;
end