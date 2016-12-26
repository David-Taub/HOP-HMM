function [startT, T, E] = genRandParams(m, n, order)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % epsilon = 0.03;
    % m x m 
    % T = eye(m) + epsilon * rand(m);
    T = rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = rand([m, ones(1,order) .* n]);
    E = bsxfun(@times, E, 1 ./ sum(E, order+1));
end

