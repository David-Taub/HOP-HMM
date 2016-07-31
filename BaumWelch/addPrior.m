
function [T1, E1] = addPrior(startT, T, E)
    [m, n] = size(E);
    E1 = [zeros(1, n); E];
    T1 = [0, startT.'; zeros(m, 1), T];
end