function [Y] = ViterbiDecode(X,T,E,startT)
% do Viterbi decoding
%
% input:
%   X:          1 x k, sequence
%   T:          m x m, transition matrix, a_ij = Prb(q_j|q_i)
%   E:          m x n, emission matrix, b_ij = Prb(o_j|q_i)
%   startT:          m x 1, prior probabilities
%
% output:
%   Y: 1 x k, sequence

	[m, n] = size(E);
	k = length(X);

	Y = zeros(1, k);
	V = zeros(m, k);     % log-prob of sequences
	ind = zeros(m, k); % records maximum state at each time stamp

	V(:, 1) = log(startT) + log(E(:, X(1)));
	for i = 2 : k
	    [V(:, i), ind(:, i)] = max(bsxfun(@plus, T', V(:, i-1)'), [], 2);
	    V(:, i) = V(:, i) + log(E(:, X(i)));
	end

	% back tracking
	[~, Y(k)] = max(V(:, k));
	for i = k-1 : -1 : 1
	    Y(i) = ind(Y(i+1), i+1);
	end
end