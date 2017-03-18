function [startT, T, E, Y, y] = genRandParamsJ(m, n, order, k)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % epsilon = 0.03;
    % m x m 
    % T = eye(m) + epsilon * rand(m);
    T = rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    
    Y = rand(m, k);
    Y = bsxfun(@times, Y, 1 ./ sum(Y, 2));

    y = rand(m, 1);
    
    % m x n
    E = rand([m, ones(1,order) .* n]);
    E = bsxfun(@times, E, 1 ./ sum(E, order+1));
end

