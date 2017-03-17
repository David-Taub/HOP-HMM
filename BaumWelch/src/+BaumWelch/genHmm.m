% generate sequence from theta (startT, T and E)
% k - length of emmision 
% T - m x m transfer matrix
% E - m x n emission matrix
function [X, Y, likelihood] = genHmm(k, startT, T, E)
    Y = zeros(1, k);
    X = zeros(1, k);
    Y(1) = find(mnrnd(1, startT));
    X(1) = find(mnrnd(1, E(Y(1), :)));
    for i = 2:k
        Y(i) = find(mnrnd(1, T(Y(i-1), :)));
        X(i) = find(mnrnd(1, E(Y(i-1), :)));
    end
end
