% Y - 1 x k
% X - 1 x k
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
function s = silhouette(X, Y, startT, T, E)
	k = length(X);
    l = X(Y==1);
	[m, n] = size(E);
	priors = statesPrior(startT, T)

	% k x n 
	EY = E(Y, :);
	% 1 x k- how well a letter is assigned to its state (smaller is better)
	% A = sum(vec2mat(X, n).' .* EY, 2).';
	A = sum(vec2mat(X, n).' .* EY, 2).' .* reshape(priors(Y), [1, k]);
	A = 1 - A;
	% 1 x k 
	nonSelfMap = (1 - vec2mat(Y, m));
	% M = (E(:, X) .* nonSelfMap);
	M = (E(:, X) .* nonSelfMap) .* repmat(priors, [1, k]);
	B = max(M, [] , 1);
	B = 1 - B;
	s = mean((B - A) ./ max([A; B]));
end


function ret = statesPrior(startT, T)
	[vecs, vals] = eig(T.');
	[~, ind] = max(abs(diag(vals)));
	ret = abs(vecs(:, ind));
	ret = ret ./ sum(ret);
end