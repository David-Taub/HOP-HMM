function [startT, T, E] = genRandParams(m, n)
    [startT, T, E] = genRandParamsExt(m, n);
    [T, E] = addPrior(startT, T, E);
end

